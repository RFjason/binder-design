# 冻结设计与筛选方法

本目录记录已完成的 IL‑7Rα minibinder 运行。正式轮已提前停止生成，当前候选数以 `results/completed.json` 为准。`run/` 是实际输入及脚本快照，其中保留最初五波计划配置；已执行 wave 0–2，wave 3–4 未执行。

## 环境及工具

复用思源一号已部署的 RFD4、BindCraft2 和 ESMFold2 三套隔离环境。生成、回折、结构审计及证据打包通过 Slurm 运行在计算节点。登录节点用于提交、文件传输和日志查询。

| 模块 | 冻结环境记录及用途 |
|---|---|
| RFD4‑Proteina | 0.0.3；Python 3.12；PyTorch 2.14.1+cu126；源码 commit `6fc44c060ef96d4c759bafeaa1eae1a66a5be67c`；NVIDIA 权重 revision `a06e2509c7033b010fc43614e9fd840303588502`。生成骨架与序列。 |
| BindCraft2 | 1.0.1；Python 3.14.7；JAX 0.11.2；使用既有 AF2 参数及 ProteinMPNN 后处理。 |
| ESMFold2 Full | 原生 `esm 3.4.1.post1`，Python 3.12、PyTorch 2.11.0+cu128；48 层 folding trunk、MSA encoder；本地完整 Full 权重。独立回折两个路线的候选。 |
| ipSAE | DunbrackLab 实现；脚本 SHA‑256 `10cf9b08c68c91e06cb28526cf2026f47a3980c9048fd3226d13e3304eaf1c27`。 |
| 结构审计 | NumPy、Biotite、tmtools 0.3.0；TM-align 骨架比较、几何审计和聚类。 |

环境预检记录及源码、权重清单位于 [provenance/environment](provenance/environment)。这些文件是复用环境的预检证据，实际每次预测的参数、驱动哈希、设备节点及软件版本保存在各候选 `summary.json`。模型权重和凭证不随结果分发。移植时需重新配置绝对路径、GPU 驱动及依赖并通过小规模回归验收。

上游工具：[RFD4](https://github.com/RosettaCommons/RFD4-Proteina)、[BindCraft2](https://github.com/PacesaLab/BindCraft2)、[ESM](https://github.com/Biohub/esm)、[ipSAE](https://github.com/DunbrackLab/IPSAE)。运行配置而非上游当前默认值定义本次实验。

## 受体与 MSA

主序列为 P16871‑1 E21–D239，219 aa。设计受体 `receptor.pdb/cif` 的 A 链按 1–219 编号，映射到 UniProt 时加 20。模型是完整野生型核心的 ESMFold2 预测；3DI3/3UP1 用于几何参考，参考结构的 T66 不替换本次野生型 I66。C 端实验 His 标签未加入本地主模型。

受体 A3M 来自公共靶标序列的 ColabFold MMseqs2 env 搜索，共 1,033 行、1,022 条不同对齐序列。复合物只给受体链 MSA；binder 链及 binder 单体不提供 MSA。FASTA、MSA、模型和冻结协议见 [run/formal](run/formal)。

- 主序列 SHA‑256：`ecdc20f8346168dbf6fd9c54ee032b1b9f1dc64a43d027bb0e06f66017de371c`
- 主 A3M SHA‑256：`7d2c71c01ea5964b735494e68e0d9d251f4d227d12cfee6d8a45158663d325d9`
- 历史 M1–G236 构建体使用独立 FASTA/A3M，只用于最终 20 条的敏感性复核，不参与主分数排序。

## 生成配置

正式轮每波 RFD4 40 个子任务、每任务 20 次采样，共 800 次；BindCraft2 60 条独立轨迹。已执行三波，即 2,400 次 RFD4 采样及 180 条 BindCraft2 轨迹。原计划五波的配置全部保留以区分“已配置”与“已执行”；执行账本见 `run/formal/state/jobs.json`。

| 表面 | UniProt 热点 | A 链编号 | 每波 RFD4 / BC2 |
|---|---|---|---:|
| 主表位 known | Y159、H211、Y212 | 139、191、192 | 400 / 30 |
| adjacent_1 | K158、Y159、K161 | 138、139、141 | 120 / 9 |
| adjacent_2 | K104、K120、T125 | 84、100、105 | 120 / 9 |
| other_1 | R86、Q89 | 66、69 | 80 / 6 |
| other_2 | E172、E175 | 152、155 | 80 / 6 |

表面选择预检要求侧链可及面积至少 24 Å²、模型 pLDDT 至少 75、可见糖链和潜在糖位点 ND2 净空至少 8 Å、末端净空至少 10 Å，并检查热点空间聚集。每个位点的具体原子选择以 `chosen_patches.json` 及逐任务配置为准。

RFD4 使用 `contig=<55-75 或 76-95>, /0, A1-219`；固定受体序列和 N/CA/C/O 主链条件；热点为指定侧链原子。参数为 `inference_experiment=quickstart`、`recycle_strategy=high/front_loaded`、`num_replicates=20`、`batch_size=1`、`num_workers=0`、`compile_lengths=null`、`metrics_tags=[]`、`prevalidate_inputs=true`。正式轮种子为 `730000 + design_batches.json 中任务索引`，以各批次保存值为准。每个 RFD4 表位组两档长度均匀分配。

BindCraft2 每任务 `max_trajectories=1`、`number_of_final_designs=1`、`design_workers=1`、`auto_multi_gpu=false`，保留失败轨迹及回折供独立统一审计。内部阈值包括 refine pLDDT 0.65、anneal/harden/mutate ipTM 0.5、final ipTM 0.6、final monomer pLDDT 0.7。主表位长度 76–95 aa；其他表面约三分之一任务 55–75 aa，其余 76–95 aa。这些工具内部阈值不代替最终统一筛选。

收集器验证恰好两条蛋白链、受体序列完全匹配、binder 长度和原子坐标有效，排除与缓存公开序列完全相同的序列；按完整 binder 序列 SHA‑256 去重，保留所有来源。新序列与 78 条小测序列再次去重。BindCraft2 内部未接受的轨迹若被收集，仍须完成相同的独立 ESMFold2 和几何筛选；因此内部 accepted 数不等于这里的通过数。

## 独立回折与界面评分

每条唯一序列均预测一次主构建体复合物及一次 binder 单体。预测器为 **ESMFold2 Full**，参数如下。

| 参数 | 数值 |
|---|---|
| seed / num_diffusion_samples | 0 / 1 |
| num_loops / num_sampling_steps | 20 / 200 |
| msa_max_depth | 1024 |
| msa_column_mask_rate / lm_dropout | 0.1 / 0.3 |
| 模板 / initial guess | 无 / 无 |

预测输出的 PAE token 顺序必须与 CIF 残基顺序核验一致，再转换成 ipSAE 输入。使用 PAE cutoff 10 Å、接触距离 cutoff 10 Å，保留两个方向，主分数为较低值 `min(A_to_B, B_to_A)`。单体与复合物 binder 的 Cα RMSD 在单体最佳叠合后计算；生成位姿一致性则在受体叠合后比较 binder，二者意义不同。

## 计算筛选阈值

所有数值和实际实现保存在 `run/formal/scripts/audit_pilot.py`，下表为本次运行内部筛选条件。

| 指标 | 通过条件 |
|---|---|
| 主 ipSAE_min | ≥0.60 |
| 单体平均 pLDDT | ≥75 |
| 单体与复合物 binder Cα RMSD | ≤4 Å |
| 受体叠合后的生成/回折 binder 位姿 RMSD | 至少一个来源 ≤5 Å |
| 来源设计的 Cα 8 Å 接触集合召回率 | 该来源 ≥0.40 |
| 回折后接触的指定热点数 | known/adjacent 至少 2，other 至少 1；小测按其实际配置 |
| 到 3DI3 可见晶体糖原子的最小距离 | ≥3.2 Å |
| 到潜在 N 糖位点 ND2 的最小距离 | ≥4.5 Å |
| 到晶体未解析末端代理位点的最小距离 | ≥3.2 Å |
| 受体/binder 重原子碰撞 | <2.2 Å 原子对为 0 |
| 半胱氨酸 | 无 Cys，或每个 SG 恰有一个 <3 Å 的 SG 配对 |
| 公开参考比对 | 不触发下述近似重复规则 |

公开参考为缓存 90 条 IL‑7Rα 序列。局部比对使用 Biotite 标准蛋白替换矩阵、gap `(-10,-1)`；identity 的分母是完整 alignment trace 行数，coverage 的分母是查询长度。触发“coverage ≥0.80 且 identity >0.70”，或“相对参考覆盖 >0.90 且 identity=1”时剔除。该初步过滤不代表已完成专利新颖性检索。binder 内部 NXS/T sequon 作为元数据记录，没有单独作为淘汰条件。

## 聚类与 Top20

所有 445 条通过者按主 ipSAE_min 降序遍历，以首个满足条件的 leader 归簇。单体 TM-align 两种长度归一化 TM-score 的较小值 ≥0.65 归入同一折叠簇；受体 5 Å 重原子接触残基集合 Jaccard ≥0.60 归入同一界面簇。这是依赖排序的 leader 聚类，不等同于两两关系的连通分量或聚类数的严格下界。

Top20 沿用冻结 `select_scale.py`，按主分数降序选择，每个已选亲本最多一条、每个折叠/界面簇组合最多一条。保留选择来源。最终的 20 条在历史 M1–G236 上按相同预测参数单独复核，分数不与主构建体混合排序。

## 证据文件与迁移

Git 中保存全部序列评分、445 条 PDB 及对应摘要、Top20 历史 PDB/摘要，以及正式轮和小测的配置、输入与脚本。Release 中的 ZIP 按 `evidence/<candidate_id>/{complex,monomer,generation,historical}` 保存完整证据，`historical` 只存在于 Top20。来源映射和原始生成文件哈希见 `provenance/passing_sources.json`。

`run/` 保留原始绝对路径用于溯源。转移运行前需重映射路径、安装兼容环境并单独准备合法获取的模型权重；直接运行旧控制器会读取已停止的状态。当前发布不启动任何新生成或重新评分。所有新增计算仍必须在计算节点运行。
