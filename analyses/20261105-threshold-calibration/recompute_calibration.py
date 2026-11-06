#!/usr/bin/env python3
"""从 Anthropic 公开数据复算本目录的阈值标定表与种子差。

用法:
    python recompute_calibration.py --data-dir <放两份公开文件的目录>

需要 cofold_predictions.parquet 与 design_summary.csv（或 .parquet）。
两份文件不进 Git，脚本会打印 SHA-256 供核对。

筛选口径与 config/scoring_protocol_v1.json 一致：
    cofolding_model = ef2full，设计靶点（target_form 为空），设计复合物化学计量，
    seed = 0（best_of_seed_selection: false）。
标签取 design_summary 的 binder_final，无标签的 120 条 Mature GDF-8 剔除。
"""
import argparse, csv, hashlib, io, statistics as st, sys
from collections import defaultdict
from pathlib import Path

# 各靶点设计时实际用的化学计量，其余按 1to1
STOICH = {n: "1to2" for n in ("15-PGDH", "VEGF-A", "Mature GDF-8", "Latent GDF-8")}
STOICH.update({"TNFa": "1to3", "Cas9": "rnp"})
THRESHOLDS = (0.30, 0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def labels(d):
    p = d / "design_summary.csv"
    if p.exists():
        rows = list(csv.DictReader(open(p, encoding="utf-8")))
    else:
        import pyarrow.parquet as pq
        p = d / "design_summary.parquet"
        rows = [dict(zip(pq.read_table(p).column_names, r))
                for r in zip(*[c.to_pylist() for c in pq.read_table(p).columns])]
    out = {}
    for r in rows:
        v = r["binder_final"]
        if v in (True, False): out[r["uuid"]] = bool(v)
        elif str(v) in ("True", "False"): out[r["uuid"]] = str(v) == "True"
    return out, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    a = ap.parse_args()
    d = Path(a.data_dir)
    import pyarrow.parquet as pq
    cof = d / "cofold_predictions.parquet"
    t = pq.read_table(cof).to_pydict()
    lab, sumpath = labels(d)

    rows = []
    for i in range(len(t["uuid"])):
        if t["cofolding_model"][i] != "ef2full": continue
        if (t["target_form"][i] or "") != "": continue
        if t["stoichiometry"][i] != STOICH.get(t["target"][i], "1to1"): continue
        u = t["uuid"][i]
        if u not in lab: continue
        rows.append((u, str(t["seed"][i]), float(t["ipsae_min"][i]), lab[u]))
    if not rows: sys.exit("筛完没有数据，检查 data-dir 与列名")

    seed0 = {u: (s, y) for u, sd, s, y in rows if sd == "0"}
    best = defaultdict(lambda: float("-inf"))
    for u, sd, s, _ in rows: best[u] = max(best[u], s)
    gap = [best[u] - seed0[u][0] for u in seed0]
    pos = [s for s, y in seed0.values() if y]
    neg = [s for s, y in seed0.values() if not y]

    print(f"cofold_predictions.parquet  sha256 {sha256(cof)}")
    print(f"{sumpath.name:26s}  sha256 {sha256(sumpath)}")
    print(f"\nseed0 n={len(seed0)}  真结合体 {len(pos)}  非结合体 {len(neg)}")
    print(f"真结合体中位 {st.median(pos):.4f}   非结合体中位 {st.median(neg):.4f}")
    print(f"每个设计的种子数 {len(rows)//len(seed0)}")
    print(f"5 种子取最大减 seed0：中位 {st.median(gap):.5f}  均值 {sum(gap)/len(gap):.5f}  "
          f">0.1 的 {sum(1 for x in gap if x > .1)} 条")
    flip = [u for u in seed0 if (seed0[u][0] >= .6) != (best[u] >= .6)]
    fp = sum(1 for u in flip if seed0[u][1])
    print(f"从 seed0 换成 5 种子取最大，0.6 线上改变过线状态 {len(flip)} 条"
          f"（真结合体 {fp}，非结合体 {len(flip)-fp}）")

    print("\n  阈值      TP    FP    FN    TN    留存率   命中率")
    for th in THRESHOLDS:
        TP = sum(1 for s, y in seed0.values() if s >= th and y)
        FP = sum(1 for s, y in seed0.values() if s >= th and not y)
        FN = len(pos) - TP; TN = len(neg) - FP
        print(f"  {th:.2f}   {TP:5d} {FP:5d} {FN:5d} {TN:5d}   "
              f"{TP/len(pos):6.1%}  {TP/(TP+FP):6.1%}" if TP + FP else "")


if __name__ == "__main__":
    main()
