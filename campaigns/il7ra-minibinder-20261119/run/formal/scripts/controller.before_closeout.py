"""Submit bounded Slurm waves and stop when 5,000 unique candidates are scored."""
from collections import Counter
import fcntl
import json
import math
import os
from pathlib import Path
import subprocess
import sys

assert os.environ.get("SLURM_JOB_ID"), "Slurm only"
ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)


def read(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else default


def dump(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


config = read(ROOT / "config/campaign.json")
state_path = ROOT / "state/controller.json"
ledger_path = ROOT / "state/jobs.json"


def submit(stage, script, arguments=(), dependency=None, array=None, partition=None, time=None, wave=None):
    command = ["sbatch", "--parsable", "--kill-on-invalid-dep=yes",
               "--job-name=il7ra5k_" + stage, "--output=logs/" + stage + "_%A_%a.out"]
    if dependency:
        command.append("--dependency=" + dependency)
    if array:
        command.append("--array=" + array)
    if partition:
        command.append("--partition=" + partition)
    if time:
        command.append("--time=" + time)
    if wave is not None:
        command.append(f"--export=ALL,SCALE_WAVE={wave}")
    command += ["scripts/" + script, *map(str, arguments)]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    job = result.stdout.strip().split(";")[0]
    assert job.isdecimal(), result.stdout
    ledger = read(ledger_path, [])
    ledger.append({"stage": stage, "wave": wave, "job": job, "command": command})
    dump(ledger_path, ledger)
    print(f"SUBMITTED {stage} wave={wave} job={job}", flush=True)
    return job


def array_for(indices, concurrency):
    assert indices
    return ",".join(map(str, indices)) + "%" + str(concurrency)


def submit_scoring(wave, collection_job=None):
    dependency = "afterok:" + collection_job if collection_job else None
    score = submit("score", "score.slurm", [f"inputs/wave_{wave}_score_%a.json"], dependency,
                   array=f"0-{config['score_shards']-1}%{config['score_concurrency']}",
                   partition=config["score_partition"], time="00:45:00", wave=wave)
    audit = submit("audit", "audit.slurm", dependency="afterany:" + score, wave=wave)
    submit("advance", "collection.slurm", ["controller"], dependency="afterany:" + audit, wave=wave)


def submit_wave(wave):
    batches = read(ROOT / "inputs/design_batches.json")
    jobs = []
    for method in ("rfd4", "bc2"):
        method_batches = [b for b in batches if b["method"] == method]
        indices = [i for i, b in enumerate(method_batches) if b["wave"] == wave]
        if indices:
            jobs.append(submit(method, "generation.slurm", [method],
                array=array_for(indices, config["generation_concurrency"][method]),
                partition=config[method + "_partition"], time="00:45:00", wave=wave))
    assert jobs
    collection = submit("collect", "collection.slurm", ["collect_scale"],
                        dependency="afterany:" + ":".join(jobs), wave=wave)
    submit_scoring(wave, collection)
    return jobs


def update_cost():
    ledger = read(ledger_path, [])
    jobs = ",".join(row["job"] for row in ledger)
    try:
        result = subprocess.run(["sacct", "-X", "-n", "-P", "-j", jobs,
                                 "-o", "JobID,State,ElapsedRaw,AllocTRES"],
                                check=True, capture_output=True, text=True, timeout=40)
        seconds = 0
        for line in result.stdout.splitlines():
            values = line.split("|")
            if len(values) >= 4 and "gres/gpu=" in values[3]:
                count = int(values[3].split("gres/gpu=")[1].split(",")[0])
                seconds += int(values[2]) * count
        dump(ROOT / "reports/slurm_cost.json", {"gpu_hours": seconds / 3600,
              "job_count": len(ledger), "sacct": result.stdout})
    except (subprocess.SubprocessError, ValueError) as exc:
        print("COST_QUERY_RETRY_LATER", repr(exc), flush=True)


def add_reserve(next_wave, deficit):
    batches = read(ROOT / "inputs/design_batches.json")
    patches = read(ROOT / "reports/chosen_patches.json")
    count = max(1, math.ceil(deficit / 16))
    cycle = ["known"] * 5 + ["adjacent_1", "adjacent_2", "adjacent_1"] + ["other_1", "other_2"]
    seed = max(b["seed"] for b in batches) + 1
    for index in range(count):
        patch = cycle[index % len(cycle)]
        low = 55 if index % 2 == 0 else 76
        template = next(b for b in batches if b["method"] == "rfd4" and
                        b["patch"] == patch and b["min_length"] == low)
        row = {**template, "id": f"rfd4_reserve_w{next_wave:02d}_{index:03d}",
               "wave": next_wave, "seed": seed + index}
        spec = read(ROOT / "inputs" / (template["id"] + ".json"))[template["id"]]
        hot = patches[patch]
        spec["conditions"]["C_HOT"]["true"] = [
            {"select": f"sel('A/*/{r['position']-20}/[{','.join(r['atom_names'])}]')"} for r in hot]
        dump(ROOT / "inputs" / (row["id"] + ".json"), {row["id"]: spec})
        batches.append(row)
    dump(ROOT / "inputs/design_batches.json", batches)
    print(f"RESERVE wave={next_wave} replicates={count*20}", flush=True)


def adapt_remaining(wave, diversity):
    if wave < 1:
        return
    groups = diversity["groups"]
    counts = {g: sum(w >= wave-1 for w in stats["parent_waves"]) for g, stats in groups.items()}
    novel = {g: sum(w >= wave-1 for w in stats["joint_cluster_first_waves"]) for g, stats in groups.items()}
    depleted = [g for g in groups if counts[g] >= 10 and novel[g] == 0]
    alternatives = [g for g in groups if novel[g] > 0]
    if not depleted or not alternatives:
        return
    destination = max(alternatives, key=lambda g: novel[g] / max(1, counts[g]))
    patches = read(ROOT / "reports/chosen_patches.json")
    destinations = [p for p in patches if p.split("_")[0] == destination]
    batches = read(ROOT / "inputs/design_batches.json")
    changed = []
    counters = Counter()
    for row in batches:
        group = row["patch"].split("_")[0]
        if row["wave"] <= wave or group not in depleted:
            continue
        counters[group] += 1
        if counters[group] % 2:
            continue
        patch = destinations[(counters[group] // 2) % len(destinations)]
        hot = patches[patch]
        spec = read(ROOT / "inputs" / (row["id"] + ".json"))
        if row["method"] == "rfd4":
            spec[row["id"]]["conditions"]["C_HOT"]["true"] = [
                {"select": f"sel('A/*/{r['position']-20}/[{','.join(r['atom_names'])}]')"} for r in hot]
        else:
            spec["targets"][0]["hotspots"] = ",".join(str(r["position"]-20) for r in hot)
        row.setdefault("original_patch", row["patch"])
        row["patch"] = patch
        row["min_hotspots"] = 1 if patch.startswith("other") else 2
        row["reallocation_after_wave"] = wave
        dump(ROOT / "inputs" / (row["id"] + ".json"), spec)
        changed.append(row["id"])
    if changed:
        dump(ROOT / "inputs/design_batches.json", batches)
        dump(ROOT / "reports" / f"reallocation_wave_{wave}.json",
             {"from": depleted, "to": destination, "tasks": changed, "new_joint_clusters": novel})


with (ROOT / "state/controller.lock").open("a") as handle:
    fcntl.flock(handle, fcntl.LOCK_EX)
    state = read(state_path, {"submitted_waves": [], "score_retries": {}})
    if len(sys.argv) > 1 and sys.argv[1] == "start":
        assert not state["submitted_waves"], "Campaign already started"
        submit_wave(0)
        state["submitted_waves"].append(0)
        state["status"] = "ACTIVE"
    else:
        wave = int(os.environ["SCALE_WAVE"])
        candidates = read(ROOT / "reports/candidates.json", [])
        missing = [r["id"] + "_" + kind for r in candidates for kind in ("complex", "monomer")
                   if not read(ROOT / "predictions" / (r["id"] + "_" + kind) / "summary.json", {}).get("ok")]
        retries = state["score_retries"].get(str(wave), 0)
        if missing:
            if retries >= 2:
                state.update(status="SCORING_NEEDS_ATTENTION", missing=missing)
                dump(state_path, state)
                raise SystemExit("Repeated prediction failures. Inspect prediction failure reports.")
            state["score_retries"][str(wave)] = retries + 1
            collection = submit("recollect", "collection.slurm", ["collect_scale"], wave=wave)
            submit_scoring(wave, collection)
        else:
            audit = read(ROOT / "reports/pilot_audit.json", {})
            diversity = read(ROOT / "reports/diversity.json", {})
            assert len(candidates) == audit.get("scored") and diversity.get("wave") == wave, "Audit failed or incomplete"
            if len(candidates) == config["new_quota"]:
                history = submit("historical", "score.slurm", ["inputs/historical_score_%a.json"],
                    array="0-3%2", partition=config["score_partition"], time="00:15:00", wave=wave)
                submit("finalize", "collection.slurm", ["finalize_scale"],
                       dependency="afterany:" + history, wave=wave)
                state["status"] = "SCORING_COMPLETE_FINAL_EXPORT_QUEUED"
            else:
                adapt_remaining(wave, diversity)
                batches = read(ROOT / "inputs/design_batches.json")
                next_wave = wave + 1
                if not any(b["wave"] == next_wave for b in batches):
                    add_reserve(next_wave, config["new_quota"] - len(candidates))
                submit_wave(next_wave)
                state["submitted_waves"].append(next_wave)
                state["status"] = "ACTIVE"
    dump(state_path, state)
    update_cost()
    print(json.dumps(state), flush=True)
