# TNFR1 扩量亲本代表复评

从本轮 97 条初筛通过序列中，对每个亲本取 v1 ipSAE_min 最高的一条。20 个代表中，2 条与先前复评的 VHH 序列及设计 ID 完全一致，复用已验证的 v2 结果。其余 18 条在 A800 作业 `64891778_0` 完成 18 个复合物和 18 个单体预测。运行时间为 2026-11-10 21:57:08 至 22:09:35 CST，用时 12 分 27 秒，退出码 `0:0`。

本批 VHH 序列保持 v1 原样，复评输入为 TNFR1 无标签 L30–T211 核心及其 MSA。预测、打分和审计沿用此前的成熟核心协议。新跑通过 **14/18**，加上复用的两条为 **16/20**。初筛通过门槛为 ipSAE_min ≥ 0.60、单体 pLDDT ≥ 70、DockQ ≥ 0.23，且单体二硫键几何通过。

| 代表序列后缀 | v1 ipSAE_min | 成熟核心 ipSAE_min | 成熟核心 DockQ | 结果 |
| --- | ---: | ---: | ---: | --- |
| RFD4 `63553d210326` | 0.702610 | 0.728803 | 0.9013 | 通过，本批最高 |
| BindCraft2 骨架的 MPNN 变体 `1da2eecfc2e5` | 0.723061 | 0.710424 | 0.8669 | 通过 |
| RFD4 骨架的 MPNN 变体 `e79afce00578` | 0.701775 | 0.000000 | 0.0094 | 未通过 |

原版评分最高的代表在成熟核心上保留界面。风险亲本 `rfd4_crd2_compact_00:tnfr1_crd2_compact_rep1_generated` 的本轮代表在 VHH 序列保持原样的情况下再次失去原设计姿态。两条通过者的 ipSAE_min 低于 0.601，分别为 `7a33037ee399` 的 0.600681 和 `7426f01b3eaa` 的 0.600885。亲本其余序列仍需根据各自的实测结果判断。

- [parent_representatives.csv](parent_representatives.csv) 合并 18 条新复评和 2 条复用结果。
- [mature_core_scores.json](mature_core_scores.json) 为作业原始报告，SHA-256 为 `c148851f28d7ff9e852684797c7ff1aedb270e76a88962f2e1fced604df26c22`。
- [provenance.json](provenance.json) 记录输入、脚本和既有快照的 SHA-256。复现时用 `prepare_parent_representatives.py <上一批暂存目录> <新暂存目录>` 重建本批输入。

历史批次另有 3 个在本轮 20 个亲本之外的亲本，其代表在[此前复评](../results/README.md)中已通过。结构去冗余、新颖性和带 His-tag 实验构建的检查仍待完成。
