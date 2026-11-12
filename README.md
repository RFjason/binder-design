# VHH 四靶点交付（FR4 修正 + CDR3 三档 + 逐条 pLDDT）

**交付时间**：2026-11-11
**平台**：AsterFire AI 工作流平台
**方法**：RFantibody 三阶段（RFdiffusion → ProteinMPNN → RF2），框架固定，只设计 CDR1/2/3

---

## 一、已落实的复核意见

| 编号 | 意见 | 落实 |
|---|---|---|
| **B1** | 补 FR4 末位 S，然后重新折叠重新打分 | ✅ 框架补至 128aa（WGQGTLVTVSS），全部重跑+重新折叠打分 |
| **B3** | CDR3 至少铺三档 | ✅ 每靶点铺 12/15/18 三档 |
| **B4** | GLP-1R 暂停 | ✅ 未提交（等组委会口径） |
| **B5** | 交付序列 + 每条 pLDDT | ✅ 本包提供 |

---

## 二、框架版本登记（供各产线统一取用）

```
文件: h-NbBCII10_128_hlt.pdb  (128 aa，FR4 末端 WGQGTLVTVSS)
MD5:    feda52b8059049971d85bb0249cb6863  (坐标文件)
SHA256: b9b6fc7b775b23a2699b21fb60f4f3237c4f20e578a0229cb878709479da2860

序列(128aa):
QVQLVESGGGLVQPGGSLRLSCAASGGSEYSYSTFSLGWFRQAPGQGLEAVAAIASMGGLTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAAVRGYFMRLPSSHNFRYWGQGTLVTVSS
序列SHA256(前32): cd82e1a37f2b55d7391ad05b8b46c3b5
```

**源头说明**：原框架取自晶体结构 3EAK（NbBCII10），该结构只解析到 127 位，末端 Ser 未解析。已按标准肽键几何（NeRF 构建，所有键长符合标准）补上第 128 位 Ser。补后：CA 链无断裂、无 clash（最近非键距离 2.255 Å 为肽键 N···O 标准距离）。

---

## 三、结果汇总

| 靶点 | 序列数 | 最高 pLDDT | 最优 CDR3 | Top 候选 | 长度 |
|---|---|---|---|---|---|
| 靶点1 EGFR | 119 | **91.73** | 15 | L02_038 | 127 aa |
| 靶点2 TNF-α | 119 | **91.93** | 15 | L02_006 | 127 aa |
| 靶点4 TfR1 | 120 | **92.33** | 15 | L02_034 | 127 aa |
| 靶点6 TNFR1 | 117 | **92.28** | 15 | L02_038 | 127 aa |
| **合计** | **475** | — | — | — | 124-130 aa |

**合规：475/475 = 100%**（判据：长度100-150aa + 2个Cys + FR1 + FR2 + FR4末端VTVSS）

### CDR3 三档分层（pLDDT）

| 靶点 | CDR3=12 | CDR3=15 | CDR3=18 |
|---|---|---|---|
| EGFR | 86.8 (max 89.8) | **89.3 (max 91.7)** | 88.7 (max 91.0) |
| TNF-α | 87.5 (max 90.2) | 88.9 (max 91.9) | **89.4 (max 91.4)** |
| TfR1 | 86.3 (max 88.3) | **89.1 (max 92.3)** | 88.6 (max 91.9) |
| TNFR1 | 87.8 (max 89.5) | **89.4 (max 92.3)** | 89.2 (max 91.6) |

**结论**：CDR3=15 在 4 个靶点中 3 个最优，TNF-α 的 18 略优。三档铺开是必要的——不同靶点最优长度不同。

---

## 四、文件说明

每个靶点目录下：

| 文件 | 内容 |
|---|---|
| `candidates.csv` | **全指标**：rank, id, design_loops, length, cdr3_len, cdr3_seq, n_cys, fr4, **plddt**, sequence |
| `candidates.fasta` | 序列（header 含 plddt / cdr3_len / loops） |
| `report.md` | 该靶点汇总报告（按 CDR3 档位统计 + Top 排行） |

---

## 五、⚠️ 未交付项与说明

| 项目 | 状态 | 原因 |
|---|---|---|
| **GLP-1R** | ❌ 未提交 | 按 B4：赛题要求激动剂，但结合 ECD 的 VHH 无激活机制；等组委会排名口径答复 |
| **TL1A / IL-7Rα**（minibinder） | ❌ 未交付 | 见下节 B2 结论 |
| **复合物预测/亲和力** | ❌ 未做 | 本包只含单体折叠 pLDDT；结合亲和力需复合物预测+对接 |
| **新颖性比对** | ❌ 未做 | 需外部已知 VHH/minibinder 库 |

### B2（minibinder）结论：骨架不可信，非偏置问题

按复核要求「先加闸门、再放宽偏置」执行后：
- 疏水核从 17.5-27.9% 提升到均值 26.6%（最大 44.1%），66.9% 达 ≥24%
- **但折叠闸门 0/8 通过**：8 个 RFdiffusion3 骨架的 ESMFold2 单体 pLDDT 全部 0.545-0.632（要求 ≥0.70）
- 同时暴露第三种组成症状：**脯氨酸 42%**（soluble 权重在低复杂度区段的偏好）

**结论**：复核人的判断被完整验证——"偏置只决定每个位置倾向填什么，不决定骨架能不能折"。RFdiffusion3 为 TL1A/IL-7Rα 生成的骨架本身不能折叠，换偏置只会换症状。minibinder 需要更换骨架生成方案（如 BindCraft2 二聚体 / RFD4-Proteina）或接受现有 WARN 级候选。

---

## 六、配置记录（B5 要求）

**RFantibody 参数**
```
input_mode: full_pipeline_structures
framework_format: hlt
framework: h-NbBCII10_128_hlt.pdb (128aa)
design_loops: scatter × 3 → "H1:10,H2:6,H3:12" / "H1:10,H2:6,H3:15" / "H1:10,H2:6,H3:18"
mpnn_loops: H1,H2,H3
num_designs: 10 (每档)
seqs_per_struct: 4
diffuser_t: 50
final_step: 1
num_recycles: 3
no_trajectory: true
scatter concurrency: 2
```

**ESMFold2 验证**
```
model_variant: fast
num_loops: 3
num_sampling_steps: 50
num_diffusion_samples: 1
output_format: pdb
```

**平台吞吐参考**（基于本轮实测）
- RFantibody 三档（10骨架×4序列×3档）：约 3-5 小时（含 GPU 排队）
- ESMFold2 折叠 ~120 条短序列：约 2-4 小时（concurrency 2）
- 单靶点全流程：约 5-6 小时

---
*由 AsterFire 平台生成*
