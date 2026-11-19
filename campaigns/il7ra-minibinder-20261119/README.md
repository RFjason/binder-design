# IL‑7Rα minibinder 设计与评分收尾结果

截至 **2026‑11‑19 05:02（北京时间）**，冻结候选池的 **2,855 条唯一序列**全部完成 ESMFold2 Full 复合物和单体评分及结构审计，**445 条通过计算初筛，20 条入选多样性 Top20**。新生成已按要求停止，控制器状态为 `CLOSEOUT_COMPLETE`。

这是 IL‑7Rα 的 **55–95 aa minibinder** 项目。仓库早期 TNFR1 VHH 的 100–150 aa 框架要求适用于 VHH 项目。

## 下载与浏览

- [全部 2,855 条评分 CSV](results/scored_candidates.csv) · [完整逐候选 JSONL](results/scored_candidates.jsonl)
- **全部 445 条通过候选**：[FASTA](passing/candidates.fasta) · [评分和结构索引 CSV](passing/candidates.csv) · [详细 JSONL](passing/candidates.jsonl) · [复合物及单体 PDB 和预测摘要](passing/structures)
- **多样性 Top20**：[FASTA](top20/candidates.fasta) · [主构建体及历史构建体评分](top20/candidates.csv) · [逐条选择及亲本来源](top20/selection.json)
- [完整原始证据 Release](https://github.com/RFjason/binder-design/releases/tag/il7ra-closeout-20261119)：445 条候选的原始生成 CIF、复合物/单体 CIF/PDB、PAE、ipSAE 输出、链与 token 映射，以及 Top20 历史构建体复核。附件校验值见 [archive_manifest.json](provenance/archive_manifest.json)。
- [设计、部署与筛选方法](METHODS.md) · [冻结运行输入及脚本](run) · [原始完成统计](results/completed.json)

## 实际设计过程

1. 冻结 P16871‑1 **E21–D239（219 aa）**序列、野生型受体模型和 1,033 行受体 A3M。设计结构 A 链 1–219 对应 UniProt 21–239；本地主模型没有 C 端实验 His 标签。3DI3 用于配体表面与可见糖链几何参考，模型保留 I66。
2. 小测完成 **78 条**去重候选的生成、回折和审计，其中 **9 条**通过，随后正式放大。正式轮以 RFD4 和 BindCraft2 两条路线探索主表位、相邻表面和其他表面，计划权重为 50% / 30% / 20%，binder 长度为 55–75 或 76–95 aa。
3. 原计划总池 5,000 条。实际执行 wave 0–2，每波 800 次 RFD4 采样、60 条 BindCraft2 轨迹；正式轮累计登记 **919 → 1,864 → 2,777 条**新唯一序列。每条轨迹可产生多个待审计序列；轨迹数和唯一序列数不是同一个统计量。
4. 将 2,777 条正式轮候选与 78 条小测候选合并，逐序列独立预测复合物和单体。保存 ESMFold2 Full 输出及 PAE/token/CIF 映射，计算双向 ipSAE、单体 pLDDT、回折/位姿一致性、糖链及末端净空、碰撞和半胱氨酸配对。
5. 对初筛通过者做公开参考序列比对，再用 TM-align 骨架聚类和受体接触集合聚类。按冻结规则选出 Top20，并在 **M1–G236 历史构建体**上复评。
6. 11 月 17 日停止后续生成，取消 advance 作业 `65615043`，保留评分、审计和收尾；最终在 11 月 19 日 05:02:05 完成导出。实际总量 2,855 条，未继续补至 5,000。

这里记录本次实际使用的构建体和阈值。赛事页面曾出现构建体版本差异，主评分与历史复评分分别保存；本结果不把任一历史页面视为当前赛事评分输入的确认。

## 结果

| 路线 | 去重并评分 | 通过计算初筛 | 通过率 |
|---|---:|---:|---:|
| RFD4 | 2,448 | 428 | 17.48% |
| BindCraft2 | 407 | 17 | 4.18% |
| 合计 | **2,855** | **445** | **15.59%** |

正式轮通过 436 条，加小测 9 条。445 条覆盖 **209 个折叠簇、18 个界面簇**；审计执行错误为 0。上述通过率描述本次目标、参数、资源配额及后处理流程。

Top20 按 `ipSAE_min` 降序遍历，每个已选亲本最多一条，每个“折叠簇、界面簇”组合最多一条，最多选 20 条。该集合含 **20 个亲本、20 个折叠簇、3 个界面簇**，均来自 RFD4。选择没有强制路线配额或界面配额。另存 [仅按分数前 20 条](top20/score_only_first20.csv) 便于比较。

主构建体 Top20 `ipSAE_min` 为 **0.830419–0.867592**，历史构建体为 **0.807653–0.872165**，20 条均在两个构建体上高于 0.6；单体 pLDDT 为 **80.63–92.88**。

| 排名 | 候选与复合物结构 | aa | 主 ipSAE_min | 历史 ipSAE_min | 单体 pLDDT | 折叠 / 界面簇 |
|---:|---|---:|---:|---:|---:|---|
| 1 | [IL7RA_s5k_fba3f5f69438d184](passing/structures/IL7RA_s5k_fba3f5f69438d184/complex/model.pdb) | 83 | 0.867592 | 0.846922 | 89.00 | 0 / 0 |
| 2 | [IL7RA_s5k_023286807b1c761e](passing/structures/IL7RA_s5k_023286807b1c761e/complex/model.pdb) | 89 | 0.866770 | 0.860150 | 89.71 | 1 / 0 |
| 3 | [IL7RA_s5k_ebf3cefdc2c6d4aa](passing/structures/IL7RA_s5k_ebf3cefdc2c6d4aa/complex/model.pdb) | 78 | 0.866547 | 0.856475 | 92.87 | 2 / 1 |
| 4 | [IL7RA_s5k_6007274ab95e86eb](passing/structures/IL7RA_s5k_6007274ab95e86eb/complex/model.pdb) | 59 | 0.865203 | 0.859430 | 84.46 | 3 / 0 |
| 5 | [IL7RA_s5k_bf8af894a3634434](passing/structures/IL7RA_s5k_bf8af894a3634434/complex/model.pdb) | 73 | 0.861151 | 0.870279 | 85.62 | 4 / 0 |
| 6 | [IL7RA_s5k_36d47dbbd91797f0](passing/structures/IL7RA_s5k_36d47dbbd91797f0/complex/model.pdb) | 79 | 0.858384 | 0.858470 | 92.88 | 5 / 0 |
| 7 | [IL7RA_pilot_e8446e8e908ae068](passing/structures/IL7RA_pilot_e8446e8e908ae068/complex/model.pdb) | 73 | 0.854147 | 0.860888 | 84.09 | 6 / 0 |
| 8 | [IL7RA_pilot_a8d46103b3f88d0f](passing/structures/IL7RA_pilot_a8d46103b3f88d0f/complex/model.pdb) | 76 | 0.851996 | 0.872165 | 82.80 | 7 / 0 |
| 9 | [IL7RA_s5k_ccfdfedf8b6a7a05](passing/structures/IL7RA_s5k_ccfdfedf8b6a7a05/complex/model.pdb) | 94 | 0.848487 | 0.838702 | 87.69 | 8 / 0 |
| 10 | [IL7RA_s5k_3c1066dd08bb77a2](passing/structures/IL7RA_s5k_3c1066dd08bb77a2/complex/model.pdb) | 71 | 0.846640 | 0.849239 | 86.77 | 9 / 0 |
| 11 | [IL7RA_s5k_9d6fa5f9f306f288](passing/structures/IL7RA_s5k_9d6fa5f9f306f288/complex/model.pdb) | 79 | 0.842931 | 0.809964 | 91.71 | 10 / 0 |
| 12 | [IL7RA_s5k_665ed6a398be852d](passing/structures/IL7RA_s5k_665ed6a398be852d/complex/model.pdb) | 68 | 0.841956 | 0.846282 | 80.63 | 11 / 2 |
| 13 | [IL7RA_s5k_9a9473c5864c0dea](passing/structures/IL7RA_s5k_9a9473c5864c0dea/complex/model.pdb) | 71 | 0.841657 | 0.841959 | 89.13 | 12 / 0 |
| 14 | [IL7RA_s5k_f38be16073dc9d23](passing/structures/IL7RA_s5k_f38be16073dc9d23/complex/model.pdb) | 72 | 0.839198 | 0.854825 | 87.93 | 13 / 0 |
| 15 | [IL7RA_s5k_810d33a1146e9f64](passing/structures/IL7RA_s5k_810d33a1146e9f64/complex/model.pdb) | 63 | 0.837972 | 0.854520 | 81.80 | 14 / 0 |
| 16 | [IL7RA_s5k_b99734f889cf7874](passing/structures/IL7RA_s5k_b99734f889cf7874/complex/model.pdb) | 90 | 0.836075 | 0.845287 | 91.31 | 15 / 0 |
| 17 | [IL7RA_s5k_d10d10468dcca4a2](passing/structures/IL7RA_s5k_d10d10468dcca4a2/complex/model.pdb) | 87 | 0.835128 | 0.807653 | 84.57 | 16 / 0 |
| 18 | [IL7RA_s5k_0d209c81c155522d](passing/structures/IL7RA_s5k_0d209c81c155522d/complex/model.pdb) | 66 | 0.832008 | 0.839875 | 82.03 | 17 / 0 |
| 19 | [IL7RA_s5k_654c60d1b011b1d1](passing/structures/IL7RA_s5k_654c60d1b011b1d1/complex/model.pdb) | 59 | 0.831091 | 0.816665 | 86.92 | 18 / 0 |
| 20 | [IL7RA_s5k_64097ef731c2808d](passing/structures/IL7RA_s5k_64097ef731c2808d/complex/model.pdb) | 84 | 0.830419 | 0.831295 | 88.38 | 19 / 0 |

## 证据范围及剩余检查

公开对照中 6 条 `binder_final=True` 的本地 ipSAE_min 为 **0.823597–0.883398**，6 条 `binder_final=False` 为 **0.852433–0.881719**，两组明显重叠，详见 [公开对照评分](results/public_controls.csv)。因此本次分数支持结构计算筛选，不能据此估计实验结合成功率或亲和力。公开对照仅用于校准，不在 2,855 条设计池内，也未作为设计序列起点。

已完成与缓存 90 条同靶标公开序列的初步相似性排查、可见晶体糖链及潜在糖位点净空审计。**完整糖链可及性、完整靶点相关专利及序列新颖性排查尚未完成；没有新候选的实验结合或功能数据。**

## 完整性检查

```bash
python tools/validate_results.py
python tools/validate_il7ra_results.py
gh release download il7ra-closeout-20261119 --repo RFjason/binder-design --dir artifacts/il7ra-closeout-20261119
python tools/validate_il7ra_results.py --archives artifacts/il7ra-closeout-20261119
```

第一个检查保留仓库已有 TNFR1 快照验证。第二个单独验证本次 IL‑7Rα 数据、筛选重算、Top20 选择、序列/结构一致性和文件哈希；带 `--archives` 时还检查 Release 证据附件。
