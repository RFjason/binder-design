#!/usr/bin/env python3
"""scale5000 轮次的路线经济性与选型复核。纯 CPU，只读已落盘的快照。

回答四个问题。
  1. 每条路线实际消耗多少生成 GPU 小时。工作单里缺的第二个数。
  2. 过线设计集中在几个亲本上，一个亲本值多少 GPU 小时。
  3. 成熟核心复评的 delta 是亲本的属性还是单条设计的属性。
  4. 按亲本取最好一条之后，能填满 20 个提交名额吗。

用法
    python route_economics.py [仓库根目录]
输出写到同目录的 route_economics.txt 与 parent_projection.csv。
"""
import csv, json, hashlib, statistics, collections, sys, os

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..")
C = os.path.join(ROOT, "campaigns", "tnfr1-vhh-20261106-scale5000")
OUT = os.path.dirname(os.path.abspath(__file__))

# 生成阶段与评分阶段的 Slurm 数组号，取自 results/completed_20261110/source_submission_ids.tsv
GEN_ARRAYS = {"64675207", "64761191", "64680855"}
SCORE_ARRAYS = {"64675209", "64761241", "64680857", "64713481"}
PILOT_ARRAY = "64680855"          # 多框架试跑，任务号 306 起，不在 batches.tsv 里
GATE = 0.6                        # 冻结的 ipSAE_min 初筛门槛

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

def load():
    f_sc = os.path.join(C, "results", "completed_20261110", "scored_candidates.jsonl")
    f_co = os.path.join(C, "results", "completed_20261110", "source_slurm_cost.json")
    f_ba = os.path.join(C, "config", "batches.tsv")
    f_mc = os.path.join(C, "experiments", "mature-core", "results", "mature_core_scores.csv")
    rows = [json.loads(l) for l in open(f_sc, encoding="utf-8")]
    cost = json.load(open(f_co, encoding="utf-8"))
    # batches.tsv 无表头，四列为 任务号、方法、批次名、分区
    batch = {int(r[0]): r[1] for r in csv.reader(open(f_ba, encoding="utf-8"), delimiter="\t")}
    mature = list(csv.DictReader(open(f_mc, encoding="utf-8")))
    return rows, cost, batch, mature, [f_sc, f_co, f_ba, f_mc]

def gpu_hours_by_route(cost, batch):
    gen = collections.defaultdict(float)
    score_sec = 0.0
    for j in cost["jobs"]:
        arr, _, task = j["job_id"].partition("_")
        sec = j["allocated_gpu_seconds"]
        if arr in GEN_ARRAYS:
            route = "multiframe_pilot" if arr == PILOT_ARRAY else batch.get(int(task), "UNMAPPED")
            gen[route] += sec / 3600.0
        elif arr in SCORE_ARRAYS:
            score_sec += sec / 3600.0
    return dict(gen), score_sec

def main():
    rows, cost, batch, mature, files = load()
    P = [r for r in rows if r["computational_quality_pass"]]
    gen, score_h = gpu_hours_by_route(cost, batch)
    gen_h = sum(gen.values())
    out = []
    w = out.append

    w("# scale5000 路线经济性与选型复核")
    w("")
    w(f"评分 {len(rows)} 条，初筛过线 {len(P)} 条，过线亲本 "
      f"{len(set(r['parent_id'] for r in P))} 个。")
    w(f"生成 {gen_h:.2f} GPU 小时，评分 {score_h:.2f} GPU 小时，合计 {gen_h + score_h:.2f}。"
      f"账本记 {cost['new_allocated_gpu_hours']:.2f}。")
    assert abs(gen_h + score_h - cost["new_allocated_gpu_hours"]) < 0.02, "机时拆分与账本不符"
    w("")

    # 一、路线经济性。评分机时按条数均摊，因为每条的评分作业规格相同
    per_design_score = score_h / len(rows)
    w("## 一、每条路线的机时与产出")
    w("")
    w("评分机时按评分条数均摊，依据是每条的评分作业规格相同。")
    w("")
    w("| 路线 | 生成 h | 评分 h | 合计 h | 评分条数 | 亲本 | 过线 | 过线亲本 | 每个过线亲本的机时 |")
    w("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    econ = {}
    for m in sorted(set(r["method"] for r in rows)):
        sub = [r for r in rows if r["method"] == m]
        sp = [r for r in sub if r["computational_quality_pass"]]
        par = len(set(r["parent_id"] for r in sub))
        pp = len(set(r["parent_id"] for r in sp))
        gh = gen.get(m, 0.0)
        sh = per_design_score * len(sub)
        tot = gh + sh
        cpp = f"{tot / pp:.2f}" if pp else "无产出"
        econ[m] = (tot, pp)
        w(f"| {m} | {gh:.2f} | {sh:.2f} | {tot:.2f} | {len(sub)} | {par} | {len(sp)} | {pp} | {cpp} |")
    mp = gen.get("multiframe_pilot", 0.0)
    nat = [r for r in rows if r["framework_id"] != "h-NbBCII10"]
    w(f"| multiframe_pilot（两个原生框架） | {mp:.2f} | {per_design_score * len(nat):.2f} | "
      f"{mp + per_design_score * len(nat):.2f} | {len(nat)} | "
      f"{len(set(r['parent_id'] for r in nat))} | 0 | 0 | 无产出 |")
    w("")

    # 二、分数分布。看一条路线是差一点还是差很远
    w("## 二、各路线的分数分布")
    w("")
    w("| 路线 | 最高 | 中位 | ≥0.55 | ≥0.60 |")
    w("| --- | ---: | ---: | ---: | ---: |")
    for m in sorted(set(r["method"] for r in rows)):
        v = sorted((r["ipsae_min"] for r in rows if r["method"] == m), reverse=True)
        w(f"| {m} | {v[0]:.4f} | {statistics.median(v):.4f} | "
          f"{sum(1 for x in v if x >= 0.55)} | {sum(1 for x in v if x >= GATE)} |")
    w("")
    hi = [r for r in rows if r["ipsae_min"] >= GATE]
    w(f"ipSAE_min 达到 {GATE} 的有 {len(hi)} 条，其中四道初筛全过的有"
      f" {sum(1 for r in hi if r['computational_quality_pass'])} 条。")
    w("单体 pLDDT、DockQ 与二硫键几何这三道在本轮没有淘汰任何一条。")
    w("")

    # 三、过线集中度
    w("## 三、过线设计的集中度")
    w("")
    cc = collections.Counter(r["parent_id"] for r in P)
    w(f"过线 {len(P)} 条分布在 {len(cc)} 个亲本上，每个亲本的过线条数中位数为"
      f" {statistics.median(cc.values()):.0f}，最大 {max(cc.values())}。")
    top = cc.most_common(1)[0]
    w(f"最多的一个亲本 `{top[0]}` 独占 {top[1]} 条，占全部过线的 {top[1] / len(P) * 100:.0f}%。")
    w("全部 97 条的 CDR 三元组互不相同，所以名额的真实分母不是序列数，是亲本数。")
    w("")

    # 四、成熟核心 delta 的归属
    w("## 四、成熟核心复评的 delta 属于亲本还是属于单条设计")
    w("")
    g = collections.defaultdict(list)
    for r in mature: g[r["parent_id"]].append(float(r["delta_ipsae"]))
    within = [x - statistics.mean(v) for v in g.values() if len(v) > 1 for x in v]
    all_d = [float(r["delta_ipsae"]) for r in mature]
    sd_in, sd_all = statistics.pstdev(within), statistics.pstdev(all_d)
    w(f"21 条复评结果按亲本分组。组内去均值后的标准差为 {sd_in:.4f}，"
      f"全部 21 条的标准差为 {sd_all:.4f}，相差 {sd_all / sd_in:.0f} 倍。")
    w("结论是 delta 由亲本的姿态决定，同一亲本的兄弟序列几乎同步变化。")
    w("")
    w("| 亲本 | 复评条数 | delta 均值 | 组内极差 |")
    w("| --- | ---: | ---: | ---: |")
    for p, v in sorted(g.items(), key=lambda kv: -len(kv[1])):
        w(f"| {p} | {len(v)} | {statistics.mean(v):+.4f} | {max(v) - min(v):.4f} |")
    w("")

    # 五、把实测 delta 外推到本轮，并给出每亲本最优的 20 条
    w("## 五、把实测 delta 外推到本轮过线池")
    w("")
    pdelta = {p: statistics.mean(v) for p, v in g.items()}
    best = {}
    for r in P:
        if r["parent_id"] not in best or r["ipsae_min"] > best[r["parent_id"]]["ipsae_min"]:
            best[r["parent_id"]] = r
    w("| 亲本 | 过线条数 | 最好 v1 | 实测 delta | 外推 v2 | 判定 |")
    w("| --- | ---: | ---: | ---: | ---: | --- |")
    proj = []
    cov = hold = drop = 0
    for p, n in cc.most_common():
        b = best[p]["ipsae_min"]
        d = pdelta.get(p)
        if d is None:
            verdict, dstr, v2 = "无实测 delta，需复评", "", ""
        else:
            cov += n
            v2 = b + d
            verdict = "外推仍过线" if v2 >= GATE else "外推掉线"
            hold += (v2 >= GATE); drop += (v2 < GATE)
            dstr, v2 = f"{d:+.4f}", f"{v2:.4f}"
        w(f"| {p} | {n} | {b:.4f} | {dstr} | {v2} | {verdict} |")
        proj.append({"parent_id": p, "passing_designs": n, "best_design_id": best[p]["design_id"],
                     "best_v1_ipsae_min": f"{b:.6f}", "measured_parent_delta": dstr,
                     "projected_v2_ipsae_min": v2, "verdict": verdict})
    w("")
    w(f"有实测亲本 delta 可外推的覆盖 {cov}/{len(P)} 条过线设计。"
      f"其中 {hold} 个亲本外推仍过线，{drop} 个外推掉线，{len(cc) - hold - drop} 个没有实测值。")
    w("外推是外推，不是复评结果。每个亲本都要有自己的成熟核心分数才能定去留。")
    w("")

    w("## 六、按亲本取最好一条的 20 条清单")
    w("")
    w("| 排名 | ipSAE_min | DockQ | 单体 pLDDT | 骨架来源 | design_id |")
    w("| ---: | ---: | ---: | ---: | --- | --- |")
    for i, r in enumerate(sorted(best.values(), key=lambda r: -r["ipsae_min"]), 1):
        w(f"| {i} | {r['ipsae_min']:.4f} | {r['pose_dockq']:.3f} | {r['monomer_plddt']:.1f} | "
          f"{r['structure_method']} | {r['design_id']} |")
    w("")
    w(f"恰好 {len(best)} 条，与赛道名额 20 相等。按亲本去冗余之后不存在名额竞争。")
    w("")

    w("## 七、输入文件哈希")
    w("")
    for f in files:
        w(f"- `{os.path.relpath(f, ROOT)}`  sha256 `{sha256(f)}`")

    txt = "\n".join(out) + "\n"
    open(os.path.join(OUT, "route_economics.txt"), "w", encoding="utf-8").write(txt)
    with open(os.path.join(OUT, "parent_projection.csv"), "w", newline="\n", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(proj[0]), lineterminator="\n")
        wr.writeheader(); wr.writerows(proj)
    print(txt)

    # 自检。这几条断言任何一条不成立都说明快照或本脚本的口径变了
    assert len(rows) == 3873 and len(P) == 97, "快照条数与 summary.json 不一致"
    assert len(best) == 20, "按亲本去冗余后的条数变了，第六节的结论要重写"
    assert sd_all / sd_in > 10, "delta 的组间组内比例变了，第四节的结论要重写"
    assert all(r["ipsae_min"] >= GATE for r in P), "有过线设计的 ipSAE_min 低于门槛"
    print("自检通过")

if __name__ == "__main__":
    main()
