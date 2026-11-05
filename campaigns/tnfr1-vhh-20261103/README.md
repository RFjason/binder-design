# TNFR1 VHH：2026-11-03 首轮计算快照

**状态：计算筛选结果；最终提交资格和实验结合均未验证。**

## 已完成的规模

| 路线 | 已统一评分 | 基础计算筛选通过 |
| --- | ---: | ---: |
| RFD4 | 520 | 4 |
| RFantibody | 357 | 0 |
| BindCraft2 | 250 | 15 |
| 合计 | 1,127 | 19 |

19条初筛候选覆盖8条父设计轨迹。序列变体数、父设计数和结构聚类数是不同统计口径；结构聚类终审尚未完成。

## 目录

```text
config/                  冻结评分协议
inputs/                  靶标序列/结构、MSA、框架、热点和生成配置
results/all_candidates.csv
results/all_candidates.jsonl
results/summary.json
shortlist/candidates.fasta
shortlist/candidates.csv
shortlist/structures/<candidate_id>/
provenance/               来源审计、作业台账、文件和附件SHA-256
scripts/hpc_snapshot/     历史运行脚本，保留原路径及作业信息
```

## 如何解释结果

- `computational_quality_pass`仅表示ipSAE≥0.6、单体pLDDT≥70、自一致性DockQ≥0.23和单体二硫键几何通过。
- `final_rank_eligible`当前全部为false。公开序列原创性、完整靶标糖链可及性和最终多样性选择仍未完成。
- pLDDT和ipSAE均为模型指标，不能转换成实验亲和力。自一致性DockQ比较设计模型与独立重预测模型。
- 部分候选的平均单体pLDDT通过，但CDR局部置信度偏低。不能只按ipSAE排序选用。
- 部分VHH含潜在N-X-S/T基序。该基序是序列审计线索，不等于已经发生糖基化。
- RFantibody8条候选的设计参考无法形成DockQ可评估界面，保存空值并标为不通过；其主预测结果仍保留。

## 冻结评分

ESMFold2 Full、seed0、1 sample、无模板、无初始猜测、靶标MSA、binder无MSA。主评分是双向ipSAE较小值。靶标主评分序列为UniProt P19438-1的1–211；设计背景为成熟胞外域30–211。

详细参数、目标/MSA哈希见`config/scoring_protocol_v1.json`。CSV方便查看，JSONL保存完整记录和来源。

作业台账保留不同时间点的历史状态。当前候选数量与筛选结果以`results/summary.json`和`results/all_candidates.jsonl`为准。

## 原始附件

所有1127条候选的原始生成结构、设计参考、复合物与单体CIF、PAE和输入/token映射随仓库Release发布。`provenance/archive_manifest.json`给出文件名、大小和SHA-256；`provenance/data_manifest.json`用于验证本目录不可变快照。

19条初筛候选的设计参考和预测结构直接跟随Git，日常协作不必下载全部附件。

## 来源与环境

靶标资料：[UniProt P19438](https://www.uniprot.org/uniprotkb/P19438/entry)、[PDB 9S3I](https://www.rcsb.org/structure/9S3I)、[AlphaFold DB](https://alphafold.ebi.ac.uk/entry/P19438)。

设计工具：[RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina)、[RFantibody](https://github.com/RosettaCommons/RFantibody)、[BindCraft2](https://github.com/PacesaLab/BindCraft2)。第三方工具、数据库和模型的授权以各自上游声明为准；本仓库不分发其模型权重。

HPC脚本是历史快照，包含原账户路径与旧作业号。适配其他环境时新建配置和新轮次，不直接执行历史恢复/取消/调度脚本。
