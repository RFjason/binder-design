# TNFR1 候选的提交合规与新颖性审计（2026-11-05）

审计对象：`campaigns/tnfr1-vhh-20261103/shortlist/candidates.csv` 的 19 条候选。
输入文件哈希见 `input_hashes.txt`。本次审计不重跑任何蛋白设计推理，只做序列层面的规则检查。

## 结论速览

| 检查项 | 结果 |
|---|---|
| 长度 100—150 aa | 19/19 通过（115 至 127 aa） |
| 仅 20 种标准氨基酸 | 19/19 通过 |
| 自带标签或融合片段 | 19/19 无 |
| 半胱氨酸成对 | 19/19 为 2 个，框架二硫键完整 |
| **FR4 完整性** | ☠️ **19/19 缺 IMGT 第 128 位** |
| **队内重复** | ☠️ **3 对**，按细则仅计一条，有效条数 19 → 16 |
| 对已发表结合物的新颖性 | 19/19 通过，CDR 最高一致性 47.3% |

## 一、FR4 缺一位（IMGT 128）

全部 19 条以 `WGQGTLVTVS` 收尾，而 IMGT 的 FR4-IMGT 为 118 至 128 位，
标准 VHH 以 S128 收尾，即 `...VTVSS`。

用 ANARCII 按 IMGT 编号核验（细则第四章（二）指定用 ANARCII）：

| 序列 | 末位 IMGT | 末位残基 |
|---|---|---|
| 参考框架 h-NbBCII10（取自 `inputs/h-NbBCII10_complete.cif`） | **127** | S |
| TNFR1_bulk1_rfd4_crd2_compact_00_8fced7a5fd5b | **127** | S |
| TNFR1_bc2_recovery_bc2_gpu_crd3_16_203250d5c341 | **127** | S |

**根因不在设计环节，在参考框架文件本身。** 该 cif 的 `entity_poly` 序列长 127，
已经以 `VTVS` 收尾。晶体结构的 C 端末位残基常因无序而未解析，
若框架是从坐标或该 cif 的序列字段取的，就会少掉 S128。

**影响与建议。** 细则第四章（一）要求「单个 VHH 可变结构域，**含完整框架和 CDR**」。
缺一个 C 端残基属于轻微截断，是否会被判不合规无法替组委会判断，
但补回成本极低：末尾追加一个 S，长度 115→116、127→128，仍在 100—150 aa 之内。
**建议补回，并在补回后重新打分**，因为打分用的结构是按 127 位序列折的。

此项不影响本仓库任何已记录的评分数值，仅影响提交形态。

## 二、队内重复 3 对

按细则第四章（四）与（二），比对区域覆盖总长度大于 90% 且一致性 100% 判为重复，
**抗体仅比较 CDR 区**，CDR 用 ANARCII 按 IMGT 编号确定。

```
TNFR1_bc2_recovery_bc2_gpu_crd3_16_de38b5b8dffd  与  ..._203250d5c341
TNFR1_bc2_recovery_bc2_gpu_crd3_16_9a7253825d68  与  ..._bf4208b6c1a3
TNFR1_bc2_recovery_bc2_gpu_crd2_08_60371bf2a410  与  ..._14c12a1d7ac8
```

三对均为同一父设计轨迹下的近缘变体，CDR 完全相同。
**队内重复仅计一条，所以 19 条的有效提交数是 16 条。**
这与 README 中「19 条初筛序列来自 8 条父设计轨迹」的记述一致，
但按赛事口径折算出的有效条数此前未记录。

## 三、新颖性：未发现问题

参考库为 `reference_library/`，从 RCSB 拉取七个赛题靶点的全部 580 个 PDB 条目，
抽取其中**非靶点**的蛋白链（长度 40 至 300 aa），即已发表的结合物，共 127 条。
比对用 mmseqs2，VHH 按细则口径只比 CDR 拼接。

最高命中仅 47.3% 一致性、52.8% 覆盖，远低于 100% 的判重阈值：

| 候选 | 最相似参考 | 一致性 | 覆盖 |
|---|---|---|---|
| TNFR1_bulk1_rfd4_crd3_long_01_ab9de9edd9d5 | 3G5Y_2（EGFR 相关） | 47.3% | 52.8% |
| TNFR1_bc2_recovery_bc2_gpu_crd3_12_083b9fc2dcc | 4UV7_3（GC1118A） | 27.2% | 59.5% |
| TNFR1_bc2_recovery_bc2_gpu_crd2_08_14c12a1d7ac8 | 6P50_2（anti-IL-7Rα） | 23.8% | 56.8% |
| TNFR1_bc2_recovery_bc2_gpu_crd2_08_60371bf2a410 | 6P50_2（anti-IL-7Rα） | 23.8% | 56.8% |

⚠️ **适用范围**：本参考库只覆盖**有晶体结构**的已发表结合物。
细则所述的官方参考库还包含该靶点已发表论文与公开专利中的序列，
因此本结论是下限，不等于通过官方新颖性审查。
`candidates.csv` 的 `novelty_status` 可由 `NOT_RUN` 更新为
`PDB_STRUCTURE_SUBSET_PASS`，不建议直接写成 `PASS`。

## 四、工具

`tools/check_submission.py`，纯 Python，无需 GPU。

```bash
python tools/check_submission.py <序列.csv 或 .fasta> {vhh|minibinder} \
  --refs analyses/20261105-compliance-novelty-audit/reference_library/*.fasta \
  --mmseqs <mmseqs 路径> --out report.csv
```

检查项：长度、字母表、标签与融合片段、疑似信号肽、VHH 框架完整性（FR1/FR2/FR4）、
半胱氨酸奇偶、同聚物连段、队内重复（抗体只比 CDR）、对参考库的新颖性。
CDR 优先用 ANARCII 按 IMGT 编号，装不上时退回正则启发式。

依赖：`anarcii`（可选）、`mmseqs2`（可选，不给 `--refs` 时不需要）。

## 五、与本仓库冻结评分协议的一致性说明

我方独立复刻了赛事入围闸门，所得配置与 `config/scoring_protocol_v1.json` 的
`predictor` 字段完全一致（seed 0、num_diffusion_samples 1、num_loops 20、
num_sampling_steps 200、msa_max_depth 1024、msa_column_mask_rate 0.1、
lm_dropout 0.3、无模板、仅靶点 MSA）。两边独立推导得到同一套参数，互为佐证。

复刻保真度已用公开数据验证：在 TNF-α 上照 Anthropic 已公开的 ef2full 设置重算 36 条，
与其公布的 `ipsae_min_ef2full` 秩相关 **0.974**，中位绝对差 0.011。

☠️ **一处口径差异，对本赛道无影响但对其他赛道有。**
Anthropic 的方法文档定义 ipSAE 的 target mask 为「其余全部蛋白残基」，
即多聚体靶点合为一个整体算一次，`ipsae_min` 取两个方向的最小值。
TNFR1 是单链靶点，怎么算都一样；但 TNF-α、TL1A 的三聚体与 TfR1 的二聚体必须按此处理，
逐条靶点链两两配对再聚合会得到不同的值。若后续开这些赛道，需在协议里写明。
