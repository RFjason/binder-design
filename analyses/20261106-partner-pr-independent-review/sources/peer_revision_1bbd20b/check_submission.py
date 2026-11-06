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
    """退路：用两个框架半胱氨酸与 FR4 定位。只能切出 CDR1 与 CDR3，切不出 CDR2。
    失败返回 None。"""
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


def cdr_key(seq, track):
    """抗体的唯一可比区域，三个 CDR 全齐才算成立。返回 (key, 口径)。

    口径为 "IMGT3" 表示三段齐全，按细则只比 CDR 区；
    "FULL" 表示编号不全，退回全长比较。
    ☠️ 不要只拿 CDR1 加 CDR3 当键。CDR2 单独变化是 MPNN 的常见输出，
    漏掉它会把近缘变体误判成一致性 100% 的重复。"""
    if track != "vhh": return seq, "FULL"
    c = cdrs(seq)
    if c and c.get("CDR1") and c.get("CDR2") and c.get("CDR3"):
        return "|".join((c["CDR1"], c["CDR2"], c["CDR3"])), "IMGT3"
    return seq, "FULL"


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
        key, _ = cdr_key(seq, track)
        if key in public:
            errs.append(f"与已公开设计重复（{public[key]}），按细则跨队去重会判给先提交方")
    return errs, warns, info


def dedup_within(items, track):
    """细则：覆盖总长度 > 90% 且一致性 = 100% 视为重复；抗体仅比较 CDR 区。
    返回 (重复对, 退回全长比较的条数)。"""
    dups = []
    keys = {}
    nfull = 0
    for sid, seq in items:
        k, mode = cdr_key(seq, track)
        if track == "vhh" and mode == "FULL": nfull += 1
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
    return uniq, nfull


def mmseqs_novelty(items, refs, mmseqs, workdir, track="minibinder"):
    """用 mmseqs 查与参考库的相似度。细则口径：覆盖 >90% 且一致性 100% 判重复。
    顺带报高相似但未达标的，那类有被认定为“轻微改造”的风险。

    返回 (status, best, note)。status 的四个取值互不混用：
      NOT_RUN  后端不可用或参考库不存在，根本没比
      ERROR    比对过程失败
      NO_HIT   比完了，没有命中
      HIT      比完了，有命中
    ☠️ 调用方必须区分 NOT_RUN 与 NO_HIT。没跑过的新颖性不能算通过。"""
    import tempfile, subprocess, shutil
    if not shutil.which(mmseqs) and not os.path.exists(mmseqs):
        return "NOT_RUN", {}, f"mmseqs 后端不可用（{mmseqs}）"
    refs = [r for r in refs if r and os.path.exists(r)]
    if not refs:
        return "NOT_RUN", {}, "参考库文件一个都不存在"
    d = tempfile.mkdtemp(dir=workdir)

    def group(pairs):
        """按可比口径分组。IMGT3 与 FULL 两种口径的序列绝不混进同一次比对。"""
        g = {"IMGT3": [], "FULL": []}
        for sid, s in pairs:
            if not s: continue
            k, mode = cdr_key(s, track)
            if k: g[mode].append((sid, k))
        return g

    qg = group(items)
    rpairs = [(f"{os.path.basename(rp)}::{sid}", s) for rp in refs for sid, s in read_seqs(rp)]
    rg = group(rpairs)
    note = []
    if track == "vhh":
        if rg["FULL"]:
            note.append(f"参考库有 {len(rg['FULL'])} 条无法按 IMGT 编号（非抗体或编号失败），"
                        f"只能走全长比对，结果另行标注")
        if qg["FULL"]:
            note.append(f"候选有 {len(qg['FULL'])} 条无法按 IMGT 编号，同样走全长比对")

    def search(qp, tp, tag):
        if not qp or not tp: return "NO_HIT", {}
        q = os.path.join(d, f"q{tag}.fa"); t = os.path.join(d, f"t{tag}.fa")
        with open(q, "w") as f:
            for sid, k in qp: f.write(f">{sid}\n{k}\n")
        with open(t, "w") as f:
            for sid, k in tp: f.write(f">{sid}\n{k}\n")
        out = os.path.join(d, f"hits{tag}.m8")
        cmd = [mmseqs, "easy-search", q, t, out, os.path.join(d, f"tmp{tag}"),
               "--format-output", "query,target,pident,qcov,tcov,evalue", "-s", "7.5",
               "--max-seqs", "500", "-e", "10000", "-v", "0",
               "--comp-bias-corr", "0", "-k", "5"]
        try: subprocess.run(cmd, capture_output=True, timeout=1800, check=True)
        except Exception as e: return f"ERROR:{type(e).__name__}", {}
        hits = {}
        for line in open(out):
            p = line.split("\t")
            if len(p) < 6: continue
            qid, tid, pid, qcov = p[0], p[1], float(p[2]), float(p[3])
            pid = pid * 100 if pid <= 1.0 else pid
            qcov = qcov * 100 if qcov <= 1.0 else qcov
            label = f"{tid}[全长口径]" if (tag == "full" and track == "vhh") else tid
            cur = hits.get(qid)
            if cur is None or (pid, qcov) > (cur[1], cur[2]): hits[qid] = (label, pid, qcov)
        return "HIT" if hits else "NO_HIT", hits

    best, errs = {}, []
    # 同口径对同口径。另外把 FULL 口径的参考也拿去比全部候选的全长，不漏真重复。
    fq = dict(qg["FULL"])
    if track == "vhh" and rg["FULL"]:
        for sid, s in items: fq.setdefault(sid, s)
    full_q = list(fq.items())
    for qp, tp, tag in (("IMGT3", "IMGT3", "cdr"), (None, "FULL", "full")):
        q = qg[qp] if qp else full_q
        st, hits = search(q, rg[tp], tag)
        if st.startswith("ERROR"): errs.append(st); continue
        for k, v in hits.items():
            cur = best.get(k)
            if cur is None or (v[1], v[2]) > (cur[1], cur[2]): best[k] = v
    shutil.rmtree(d, ignore_errors=True)
    if errs and not best:
        return "ERROR", {}, "；".join(errs + note)
    if errs: note.append("部分比对失败：" + "；".join(errs))
    return ("HIT" if best else "NO_HIT"), best, "；".join(note)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seqfile"); ap.add_argument("track", choices=list(RULES))
    ap.add_argument("--public"); ap.add_argument("--out")
    ap.add_argument("--refs", nargs="*", default=[], help="新颖性参考库 FASTA，可多个")
    ap.add_argument("--mmseqs", default="mmseqs")
    ap.add_argument("--no-anarcii", action="store_true", help="跳过 ANARCII，直接用正则切 CDR")
    a = ap.parse_args()

    items = read_seqs(a.seqfile)
    pubitems = read_seqs(a.public) if (a.public and os.path.exists(a.public)) else []
    refitems = [p for rp in a.refs if rp and os.path.exists(rp) for p in read_seqs(rp)]
    print(f"读入 {len(items)} 条，赛道 {a.track}，长度限 {RULES[a.track][0]}-{RULES[a.track][1]} aa")
    if pubitems: print(f"公开设计语料 {len(pubitems)} 条，用于查跨队重复")
    if a.track == "vhh" and not a.no_anarcii:
        # ☠️ 候选、公开语料、参考库三方必须一起编号。
        # 只编号一方会让比对两端的可比区域定义不一致，比出来的一致性没有意义。
        allseq = [s for _, s in items] + [s for _, s in pubitems] + [s for _, s in refitems]
        anarcii_number_all(sorted({s for s in allseq if 80 <= len(s) <= 200}))
        if _CDR_CACHE: print(f"ANARCII/IMGT 编号成功 {len(_CDR_CACHE)} 条")
        if refitems:
            nref_ok = sum(1 for _, s in refitems if cdr_key(s, "vhh")[1] == "IMGT3")
            print(f"参考库 {len(refitems)} 条，其中 {nref_ok} 条编号成功、"
                  f"{len(refitems)-nref_ok} 条无法编号")
    public = {}
    for sid, s in pubitems:
        public[s] = sid
        k, mode = cdr_key(s, a.track)
        if mode == "IMGT3": public[k] = sid

    rows, nfail, nwarn = [], 0, 0
    for sid, seq in items:
        errs, warns, info = check(sid, seq, a.track, public)
        if errs: nfail += 1
        elif warns: nwarn += 1
        rows.append(dict(design_id=sid, 结果="FAIL" if errs else ("WARN" if warns else "PASS"),
                         长度=info.get("len"), Cys=info.get("Cys"),
                         疏水比例=info.get("疏水比例"), CDR3长度=info.get("CDR3_len", ""),
                         问题="; ".join(errs), 提示="; ".join(warns),
                         新颖性="NOT_RUN", 最相似参考="", **{"一致性%": "", "覆盖%": ""}))
    for r in rows:
        if r["结果"] != "PASS":
            print(f"  [{r['结果']}] {r['design_id']} len={r['长度']}  {r['问题']}{' | ' + r['提示'] if r['提示'] else ''}")

    nov_broken = False
    if a.refs:
        status, nov, note = mmseqs_novelty(
            items, a.refs, a.mmseqs,
            os.path.dirname(os.path.abspath(a.seqfile)) or ".", a.track)
        if note: print(f"\n新颖性说明：{note}")
        if status in ("NOT_RUN", "ERROR"):
            # ☠️ 请求了新颖性检查但后端没跑成，不能当成通过。
            nov_broken = True
            print(f"\n☠️ 新颖性比对状态 {status}，{note or '原因未知'}")
            print("   基础序列检查的结论仍然有效，但新颖性一项为未完成，不得写成 PASS。")
        else:
            for r in rows: r["新颖性"] = status
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

    dups, nfull = dedup_within(items, a.track)
    if a.track == "vhh" and nfull:
        print(f"\n⚠️ 有 {nfull} 条无法按 IMGT 切出完整三段 CDR，这些条按全长比较去重")
    if dups:
        print(f"\n☠️ 内部重复 {len(dups)} 对（细则：覆盖 >90% 且一致性 100%，抗体只比 CDR）")
        for b, keep in dups[:10]: print(f"   {b} 与 {keep} 重复")

    print(f"\n合计 PASS {len(rows)-nfail-nwarn}  WARN {nwarn}  FAIL {nfail}  内部重复 {len(dups)} 对")
    if nov_broken: print("新颖性 未完成")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print(f"逐条结果 -> {a.out}")
    return 1 if (nfail or dups or nov_broken) else 0


def selftest():
    """回归用例。全部用注入的编号结果，不依赖 ANARCII 与 mmseqs。
    三个案例分别对应 2026-11-06 复核查出的三个真实缺陷。"""
    base = "QVQLVESGGGLVQPGGSLRLSCAAS{c1}WGWFRQAPGQGLEAVA{c2}YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC{c3}WGQGTLVTVSS"
    def mk(c1, c2, c3):
        s = base.format(c1=c1, c2=c2, c3=c3)
        _CDR_CACHE[s] = {"CDR1": c1, "CDR2": c2, "CDR3": c3}
        return s
    # 1. 只有 CDR2 不同，不是重复。旧键只拼 CDR1 与 CDR3，会误判。
    a1 = mk("GDISNITY", "IDLRNNKT", "AAGSSGAFDY")
    a2 = mk("GDISNITY", "INLRNNET", "AAGSSGAFDY")
    dups, _ = dedup_within([("a1", a1), ("a2", a2)], "vhh")
    assert dups == [], f"只差 CDR2 被误判成重复：{dups}"
    # 2. 三段全同才算重复
    b1 = mk("GDISNITY", "IDLRNNKT", "AAGSSGAFDY") + ""
    b2 = base.format(c1="GDISNITY", c2="IDLRNNKT", c3="AAGSSGAFDY") + "A"
    _CDR_CACHE[b2] = {"CDR1": "GDISNITY", "CDR2": "IDLRNNKT", "CDR3": "AAGSSGAFDY"}
    dups, _ = dedup_within([("b1", b1), ("b2", b2)], "vhh")
    assert len(dups) == 1, f"三段全同应判重复：{dups}"
    # 3. 编号不全退回全长，不得只拿两段拼键
    c1s = base.format(c1="GDISNITY", c2="IDLRNNKT", c3="AAGSSGAFDY").replace("IDLRNNKT", "IDLRNNKK")
    _CDR_CACHE[c1s] = {"CDR1": "GDISNITY", "CDR3": "AAGSSGAFDY"}
    k, mode = cdr_key(c1s, "vhh")
    assert mode == "FULL" and k == c1s, "缺 CDR2 时必须退回全长比较"
    # 4. 后端缺失必须报 NOT_RUN，不能静默当成没命中
    st, best, note = mmseqs_novelty([("a1", a1)], [__file__], "/nonexistent/mmseqs", ".", "vhh")
    assert st == "NOT_RUN" and best == {} and note, f"后端缺失应为 NOT_RUN：{st}"
    print("selftest OK：4/4")


if __name__ == "__main__":
    if "--selftest" in sys.argv: selftest(); sys.exit(0)
    sys.exit(main())
