#!/usr/bin/env python3
"""提交合规校验。按「长江杯」细则第二节逐条检查，生成的序列先过这一关再谈打分。

用法:
    python check_submission.py <序列文件> <赛道>
    序列文件支持 FASTA 或 CSV（需有 id/design_id 与 seq/sequence 列）
    赛道取值 minibinder 或 vhh

可选:
    --public <fasta>   与已公开设计比对，命中则有被别队抢先提交的风险
    --out <csv>        把逐条结果写成表格

退出码 0 表示全部通过，1 表示有 FAIL。
"""
import sys, csv, re, argparse, os
from collections import Counter

AA = set("ACDEFGHIKLMNPQRSTVWY")
RULES = {"minibinder": (40, 100), "vhh": (100, 150)}

# 常见标签与融合片段，细则明确禁止自带
TAGS = {
    "His6": "HHHHHH", "His8": "HHHHHHHH", "FLAG": "DYKDDDDK", "Myc": "EQKLISEEDL",
    "HA": "YPYDVPDYA", "StrepII": "WSHPQFEK", "V5": "GKPIPNPLLGLDST",
    "SUMO尾": "QIGG", "TEV位点": "ENLYFQ", "凝血酶位点": "LVPRGS",
    "Fc铰链": "EPKSCDKTHTCPPCP", "Avi标签": "GLNDIFEAQKIEWHE",
}
# VHH 框架标志
VHH_FR1 = re.compile(r"^[QEDA][VLIM][QKE]L[VQLKE]")
VHH_FR2_W = re.compile(r"W[FYVG]RQ")            # FR2 起始，VHH 常见 WFRQ/WYRQ
VHH_FR4 = re.compile(r"WG[QKRAEP]G[TAS][QLPMTRSVAE][VLM]TVS[SAP]")


def read_seqs(path):
    out = []
    if path.lower().endswith((".fa", ".fasta", ".faa")):
        sid, buf = None, []
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                if sid: out.append((sid, "".join(buf)))
                sid, buf = line[1:].split()[0], []
            elif line: buf.append(line)
        if sid: out.append((sid, "".join(buf)))
    else:
        for r in csv.DictReader(open(path)):
            k = next((c for c in ("design_id", "id", "name") if c in r), None)
            s = next((c for c in ("seq", "sequence", "binder_seq") if c in r), None)
            if k and s: out.append((r[k], r[s].strip().upper()))
    return out


_CDR_CACHE = {}
IMGT = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}


def anarcii_number_all(seqs):
    """官方口径：ANARCII 按 IMGT 编号切 CDR（细则第四章（二）指定）。
    模型加载一次要一分多钟，所以**一次把全部序列批量编号**，之后只查表。"""
    if not seqs: return
    try:
        from anarcii import Anarcii
        m = Anarcii(seq_type="antibody", batch_size=max(8, min(64, len(seqs))),
                    cpu=True, verbose=False)
        res = m.number([(f"s{i}", s) for i, s in enumerate(seqs)])
    except Exception as e:
        print(f"  （ANARCII 不可用，CDR 退回正则启发式：{type(e).__name__}）")
        return
    recs = list(res.values()) if hasattr(res, "values") else list(res)
    for s, rec in zip(seqs, recs):
        try:
            num = rec["numbering"] if hasattr(rec, "keys") else rec
            out = {}
            for name, (lo, hi) in IMGT.items():
                out[name] = "".join(aa for pos, aa in num
                                    if lo <= (pos[0] if isinstance(pos, (tuple, list)) else pos) <= hi
                                    and aa != "-")
            if out.get("CDR3"): _CDR_CACHE[s] = out
        except Exception:
            pass


def cdrs(seq):
    """切 CDR。优先查 ANARCII 的批量结果，没有则退回正则启发式。"""
    if seq in _CDR_CACHE: return _CDR_CACHE[seq]
    return _cdrs_regex(seq)


def _cdrs_regex(seq):
    """退路：用两个框架半胱氨酸与 FR4 定位。失败返回 None。"""
    cys = [i for i, c in enumerate(seq) if c == "C"]
    m4 = VHH_FR4.search(seq)
    if len(cys) < 2 or not m4: return None
    c1, c2 = cys[0], max(c for c in cys if c < m4.start())
    if c2 <= c1: return None
    m2 = VHH_FR2_W.search(seq, c1)
    if not m2: return None
    cdr1 = seq[c1 + 1:m2.start()]
    cdr3 = seq[c2 + 1:m4.start()]
    return {"CDR1": cdr1, "CDR3": cdr3}


def check(sid, seq, track, public):
    errs, warns, info = [], [], {}
    lo, hi = RULES[track]
    info["len"] = len(seq)

    if not seq:
        return ["序列为空"], warns, info
    bad = sorted(set(seq) - AA)
    if bad: errs.append(f"含非标准氨基酸 {''.join(bad)}")
    if not (lo <= len(seq) <= hi):
        errs.append(f"长度 {len(seq)} 超出 {track} 赛道的 {lo}-{hi} aa")

    for name, motif in TAGS.items():
        if motif in seq: errs.append(f"疑似自带{name}（{motif}）")
    # 末端连续 His 也算标签
    if re.search(r"H{5,}$", seq) or re.match(r"^M?H{5,}", seq):
        errs.append("末端连续组氨酸，疑似 His 标签")
    # 疑似信号肽：前 25 位高度疏水且以 M 起始
    head = seq[:25]
    if seq[0] == "M" and sum(c in "AVLIMFWC" for c in head) / len(head) > 0.6:
        warns.append("N 端 25 位疏水比例偏高，确认不是信号肽")

    ncys = seq.count("C")
    info["Cys"] = ncys
    if track == "vhh":
        if ncys < 2: errs.append(f"VHH 只有 {ncys} 个半胱氨酸，框架二硫键不完整")
        elif ncys % 2: warns.append(f"半胱氨酸 {ncys} 个为奇数，有游离巯基风险")
        miss = []
        if not VHH_FR1.search(seq): miss.append("FR1 起始")
        if not VHH_FR2_W.search(seq): miss.append("FR2 的 W(F/Y)RQ")
        if not VHH_FR4.search(seq): miss.append("FR4 的 WGxGTxVTVSS")
        if miss: errs.append("VHH 框架不完整，缺 " + "、".join(miss))
        c = cdrs(seq)
        if c:
            info["CDR3"] = c["CDR3"]
            info["CDR3_len"] = len(c["CDR3"])
            if not 5 <= len(c["CDR3"]) <= 30:
                warns.append(f"CDR3 长度 {len(c['CDR3'])} 异常")
    else:
        if ncys % 2: warns.append(f"半胱氨酸 {ncys} 个为奇数，游离巯基会影响表达与均一性")
        if ncys > 2: warns.append(f"{ncys} 个半胱氨酸偏多，minibinder 通常无二硫键")

    m = re.search(r"(.)\1{4,}", seq)
    if m: warns.append(f"同聚物连段 {m.group(0)}，序列设计质量存疑")
    hyd = sum(c in "AVLIMFWY" for c in seq) / len(seq)
    info["疏水比例"] = round(hyd, 3)
    if hyd > 0.5: warns.append(f"疏水比例 {hyd:.2f} 偏高，聚集风险")
    low = Counter(seq).most_common(1)[0]
    if low[1] / len(seq) > 0.3: warns.append(f"单一残基 {low[0]} 占 {low[1]/len(seq):.0%}")

    if public:
        key = info.get("CDR3") if track == "vhh" and info.get("CDR3") else seq
        if key in public:
            errs.append(f"与已公开设计重复（{public[key]}），按细则跨队去重会判给先提交方")
    return errs, warns, info


def dedup_within(items, track):
    """细则：覆盖总长度 > 90% 且一致性 = 100% 视为重复；抗体仅比较 CDR 区。"""
    dups = []
    keys = {}
    for sid, seq in items:
        k = seq
        if track == "vhh":
            c = cdrs(seq)
            if c: k = c["CDR1"] + "|" + c["CDR3"]
        if k in keys: dups.append((sid, keys[k]))
        else: keys[k] = sid
    # 子串型重复（覆盖 >90% 且完全一致）
    for i, (sid, s) in enumerate(items):
        for sid2, s2 in items[i + 1:]:
            if s == s2: continue
            short, long_ = (s, s2) if len(s) < len(s2) else (s2, s)
            if short in long_ and len(short) / len(long_) > 0.9:
                dups.append((sid2, sid))
    seen, uniq = set(), []
    for b, keep in dups:
        k = tuple(sorted((b, keep)))
        if k not in seen: seen.add(k); uniq.append((b, keep))
    return uniq


def mmseqs_novelty(items, refs, mmseqs, workdir, track="minibinder"):
    """用 mmseqs 查与参考库的相似度。细则口径：覆盖 >90% 且一致性 100% 判重复。
    顺带报高相似但未达标的，那类有被认定为“轻微改造”的风险。"""
    import tempfile, subprocess, shutil
    if not shutil.which(mmseqs) and not os.path.exists(mmseqs): return {}
    refs = [r for r in refs if r and os.path.exists(r)]
    if not refs: return {}
    d = tempfile.mkdtemp(dir=workdir)
    q = os.path.join(d, "q.fa"); t = os.path.join(d, "t.fa")

    def key(s):
        """细则：抗体仅比较 CDR 区。VHH 框架天生保守，拿全长比会把每条都报成高相似。"""
        if track != "vhh": return s
        c = cdrs(s)
        return (c["CDR1"] + c.get("CDR2", "") + c["CDR3"]) if c else s

    with open(q, "w") as f:
        for sid, s in items:
            k = key(s)
            if k: f.write(f">{sid}\n{k}\n")
    with open(t, "w") as f:
        for rp in refs:
            for sid, s in read_seqs(rp):
                k = key(s)
                if k: f.write(f">{os.path.basename(rp)}::{sid}\n{k}\n")
    out = os.path.join(d, "hits.m8")
    cmd = [mmseqs, "easy-search", q, t, out, os.path.join(d, "tmp"),
           "--format-output", "query,target,pident,qcov,tcov,evalue", "-s", "7.5",
           "--max-seqs", "500", "-e", "10000", "-v", "0",
           "--comp-bias-corr", "0", "-k", "5"]
    try: subprocess.run(cmd, capture_output=True, timeout=1800, check=True)
    except Exception as e:
        print(f"  （mmseqs 比对跳过：{type(e).__name__}）"); return {}
    best = {}
    for line in open(out):
        p = line.split("\t")
        if len(p) < 6: continue
        qid, tid, pid, qcov = p[0], p[1], float(p[2]), float(p[3])
        pid = pid * 100 if pid <= 1.0 else pid
        qcov = qcov * 100 if qcov <= 1.0 else qcov
        cur = best.get(qid)
        if cur is None or (pid, qcov) > (cur[1], cur[2]): best[qid] = (tid, pid, qcov)
    shutil.rmtree(d, ignore_errors=True)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seqfile"); ap.add_argument("track", choices=list(RULES))
    ap.add_argument("--public"); ap.add_argument("--out")
    ap.add_argument("--refs", nargs="*", default=[], help="新颖性参考库 FASTA，可多个")
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--no-anarcii", action="store_true", help="跳过 ANARCII，直接用正则切 CDR")
    a = ap.parse_args()

    items = read_seqs(a.seqfile)
    public = {}
    if a.public and os.path.exists(a.public):
        for sid, s in read_seqs(a.public):
            public[s] = sid
            c = cdrs(s)
            if c: public[c["CDR3"]] = sid
    print(f"读入 {len(items)} 条，赛道 {a.track}，长度限 {RULES[a.track][0]}-{RULES[a.track][1]} aa")
    if public: print(f"公开设计语料 {len(public)} 条，用于查跨队重复")
    if a.track == "vhh" and not a.no_anarcii:
        allseq = [s for _, s in items] + [s for s in public if public]
        anarcii_number_all(sorted({s for s in allseq if 80 <= len(s) <= 200}))
        if _CDR_CACHE: print(f"ANARCII/IMGT 编号成功 {len(_CDR_CACHE)} 条")

    rows, nfail, nwarn = [], 0, 0
    for sid, seq in items:
        errs, warns, info = check(sid, seq, a.track, public)
        if errs: nfail += 1
        elif warns: nwarn += 1
        rows.append(dict(design_id=sid, 结果="FAIL" if errs else ("WARN" if warns else "PASS"),
                         长度=info.get("len"), Cys=info.get("Cys"),
                         疏水比例=info.get("疏水比例"), CDR3长度=info.get("CDR3_len", ""),
                         问题="; ".join(errs), 提示="; ".join(warns),
                         最相似参考="", **{"一致性%": "", "覆盖%": ""}))
    for r in rows:
        if r["结果"] != "PASS":
            print(f"  [{r['结果']}] {r['design_id']} len={r['长度']}  {r['问题']}{' | ' + r['提示'] if r['提示'] else ''}")

    if a.refs:
        nov = mmseqs_novelty(items, a.refs, a.mmseqs, os.path.dirname(os.path.abspath(a.seqfile)) or ".", a.track)
        if nov:
            hard = [(q, v) for q, v in nov.items() if v[1] >= 99.5 and v[2] > 90]
            soft = [(q, v) for q, v in nov.items() if 70 <= v[1] < 99.5]
            print(f"\n新颖性比对：{len(nov)} 条有命中")
            if hard:
                print(f"  ☠️ 判重复 {len(hard)} 条（覆盖 >90% 且一致性 100%）")
                for q, v in hard[:10]: print(f"     {q} -> {v[0]} 一致性 {v[1]:.1f}% 覆盖 {v[2]:.1f}%")
            if soft:
                print(f"  ⚠️ 高相似 {len(soft)} 条，有被判“轻微改造”的风险")
                for q, v in sorted(soft, key=lambda x: -x[1][1])[:10]:
                    print(f"     {q} -> {v[0]} 一致性 {v[1]:.1f}% 覆盖 {v[2]:.1f}%")
            for r in rows:
                v = nov.get(r["design_id"])
                if v:
                    r["最相似参考"] = v[0]; r["一致性%"] = round(v[1], 1); r["覆盖%"] = round(v[2], 1)
                    if v[1] >= 99.5 and v[2] > 90:
                        if r["结果"] == "WARN": nwarn -= 1
                        if r["结果"] != "FAIL": nfail += 1
                        r["结果"] = "FAIL"
                        r["问题"] = (r["问题"] + "; " if r["问题"] else "") + f"与参考库重复（{v[0]}）"

    dups = dedup_within(items, a.track)
    if dups:
        print(f"\n☠️ 内部重复 {len(dups)} 对（细则：覆盖 >90% 且一致性 100%，抗体只比 CDR）")
        for b, keep in dups[:10]: print(f"   {b} 与 {keep} 重复")

    print(f"\n合计 PASS {len(rows)-nfail-nwarn}  WARN {nwarn}  FAIL {nfail}  内部重复 {len(dups)} 对")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print(f"逐条结果 -> {a.out}")
    return 1 if nfail or dups else 0


if __name__ == "__main__":
    sys.exit(main())
