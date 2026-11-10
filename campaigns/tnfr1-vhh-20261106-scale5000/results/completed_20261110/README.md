# TNFR1 5,000 条评分完成

2026-11-10 09:01:45 CST，Slurm 控制器记录项目达到 **5,000 条评分和审计**。其中历史基线 1,127 条，本轮新增 3,873 条。本轮有 97 条通过计算初筛，来自 20 个登记的亲本轨迹或骨架，连同历史基线累计通过 116 条。累计分配 149.23 GPU 小时，其中本轮 75.03 GPU 小时。

| 本轮新增的骨架来源 | 完成评分 |
| --- | ---: |
| RFD4 | 3,083 |
| BindCraft2 | 726 |
| RFantibody | 64 |

骨架来源包含其后续 ProteinMPNN 变体。框架分布为 h-NbBCII10 3,782 条、native_3P0G 40 条和 native_1ZVH 51 条。最高新增 ipSAE_min 为 0.723061，候选为 `TNFR1_s5k_mpnn_parent03_t04_1da2eecfc2e5`。

这套结果采用冻结 v1 输入 TNFR1 M1–T211。另行开展的 [L30–T211 成熟核心复评](../../experiments/mature-core/results/README.md) 覆盖早期的 21 条候选，其中 15 条通过复评。全量初筛名单尚待按成熟构建复评和最终筛选，最终可提交数为 0。

| 文件 | 内容 |
| --- | --- |
| [scored_candidates.csv](scored_candidates.csv) | 本轮新增 3,873 条的序列、分数、框架和亲本 |
| [scored_candidates.jsonl](scored_candidates.jsonl) | 完整逐条审计字段和源报告引用 |
| [summary.json](summary.json) | 按路线和框架统计、324 份源报告及结果文件的 SHA-256 |
| [source_campaign_finished.json](source_campaign_finished.json) | Slurm 控制器完成记录 |
| [source_slurm_cost.json](source_slurm_cost.json) | Slurm GPU 分配机时账本 |

重建快照时，使用 `scripts/export_scored_snapshot_local.py` 读取 324 份落盘审计报告和完成记录。导出时核对长度、标准氨基酸、序列 SHA-256、与基线及本轮的精确去重，以及通过标志与计算门槛是否一致。
