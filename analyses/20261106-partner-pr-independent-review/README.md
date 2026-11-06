# Binder项目同伴PR独立检查（2026-11-06）

检查对象：EthonSun的[PR #1](https://github.com/RFjason/binder-design/pull/1)（`20313414a2d98853cf66708655edda883d0f02e2`）及[PR #2](https://github.com/RFjason/binder-design/pull/2)（`7aad964e0f348fec5f8fc06735ae8a8436695e6b`）。两份PR均未合并，本次没有发布评论、修改同伴分支或变更生产作业。

结论：FR4末端补位和0.6阈值统计有证据支持；去重工具存在可复现的错误，算力集中到BC2/CRD3的建议超出了现有数据所能支持的范围。建议修改后再合并，暂不将该工具接成正式入库闸门。

## 1. [P1] 队内去重漏掉CDR2，错误报告三对重复

位置：PR #1 `tools/check_submission.py:168`。键为`CDR1 + "|" + CDR3`，即使ANARCII已给出CDR2也不使用。

在Slurm CPU作业`64679211`上使用ANARCII 2.0.8、antibody模式、IMGT编号、2 CPU独立重跑19条原序列及19条补尾变体，作业完成`ExitCode=0:0`，耗时68秒。原序列有**19种不同的完整三CDR组合**。

|PR声称的重复对（ID尾部）|第一条CDR2|第二条CDR2|
|---|---|---|
|`de38b5b8dffd` / `203250d5c341`|IDLRNNKT|INLRNNET|
|`9a7253825d68` / `bf4208b6c1a3`|IDLRNNTT|IDMRNNTT|
|`60371bf2a410` / `14c12a1d7ac8`|IDLEKKET|IDLVKNET|

三对的差异都落在独立编号后的CDR2中。“三CDR完全相同、19条因此减为16条”的结论不成立。对全部19条三CDR拼接序列进一步检查，也未发现覆盖较长序列>90%的完全一致连续片段。这不是对官方参考库或官方比对器的最终提交认证。

建议：队内去重、公库去重及部分覆盖检查统一使用同一套完整IMGT三CDR表示，记录三个区域的边界、编号成功状态及实际比对覆盖。修订PR #1/#2的重复数量及相关文字。

证据：`results/independent_numbering.json`、`results/anarcii_raw.json`、`results/cdr_local_exact_check.json`、`results/local_checks.json`。`scripts/recompute.py`还提供只改变CDR2的最小回归案例，原函数仍错误判重复。

## 2. [P1] 新颖性工具能在未执行比对时返回PASS

位置：PR #1 `tools/check_submission.py:189`及`:268–270`。

提供合法完整VHH与实际存在的`--refs`文件，故意指定不存在的MMseqs程序。原函数返回空字典，主程序没有区分“没有命中”与“根本没跑”，最终输出`PASS 1 / FAIL 0`、退出码0。比对异常也会返回同一空字典。

建议：显式区分`NOT_RUN`、`ERROR`、`RUN_COMPLETED_NO_HIT`和`HIT`；用户请求新颖性检查而后端不可用时不能返回完成审查意义上的PASS。允许保存基础序列检查结果，但不能晋级为新颖性通过。

证据：`novelty_backend_repro.fasta`、`novelty_backend_repro.csv`、`novelty_backend_repro_stdout.txt`及`scripts/recompute.py`。

## 3. [P1] 新颖性比对的查询与参考CDR处理不对称

位置：PR #1 `tools/check_submission.py:248–251`、`:196–208`。

ANARCII批量编号只覆盖候选和`--public`，没有包括`--refs`参考库。正常候选因此可得到三CDR拼接；参考序列若不恰好与已缓存序列相同，就走只抽取CDR1/CDR3的正则退路，失败则使用全长蛋白。两端不再具有一致的比对区域定义。

这会影响一致性与覆盖率，当前CSV中的低相似度不足以直接升级为`PDB_STRUCTURE_SUBSET_PASS`。参考库收集本身可保留，但应先给候选和抗体参考链一致编号，显式标记无法编号/非抗体参考，再重跑比对并保存命令、工具版本、完整命中表。论文/专利库仍未覆盖。

## 4. [P1] 将全部预算集中到CRD3的依据混入路线和亲本变体效应

位置：PR #2 `行动建议.md:99–116`。

`stratum=crd2/crd3`混合BC2和RFantibody；`crd2_compact/long`、`crd3_compact/long`仅对应RFD4。它们不是同一路线中可直接排序的六个独立表位/长度组。因此“compact/long不如基础层”不能从该汇总表推出。

更关键的是BC2序列共250条，但来自25个有保留输出的亲本轨迹；15条过线来自4个亲本。按表位拆分：

|BC2表位|有输出序列|有输出亲本|过线序列|过线亲本|
|---|---:|---:|---:|---:|
|CRD2|110|11|4|2|
|CRD3|140|14|11|2|

CRD3的11条过线没有代表11次独立成功。当前独立过线亲本数为2比2，不能用11比4的序列数来证明CRD3应独占预算。RFD4又呈现CRD2为3/260、CRD3为1/260的相反方向。

BC2值得提高优先级，但应按方法×表位分层，以所有启动轨迹为分母，包含无输出轨迹，并比较独立过线亲本/GPU小时和结构簇数量。“1127条序列预算→68条过线”不是等GPU成本预测。

证据：`results/local_checks.json`的`route`和`patch`；复算来源为冻结`all_candidates.jsonl`。

## 5. [P2] 当前过线差异不能归因于DockQ闸门

位置：PR #2 `行动建议.md:43–46`及后续机制解释。

|路线|ipSAE≥0.6|完整初筛通过|在ipSAE过线者中因DockQ失败|
|---|---:|---:|---:|
|RFD4|4|4|0|
|RFantibody|0|0|0|
|BC2|15|15|0|

DockQ分布差异是实际可观察的关联，但没有额外淘汰任何一条ipSAE过线者。“关键不在ipSAE、在DockQ”不是这批数据的门槛归因结果。按路线分组后空状态/ASSESSED中位数相近，也不能单靠分布检验证明几何实现没有问题。

此外BC2设计使用AF2-Multimer权重，而独立评分用ESMFold2；“优化过结构预测器”不等于保证被另一个预测器复现。算法或初始猜测的因果解释应保留为待检验假说。低自洽DockQ也不能替代逐例结构检查，更不是实验结合失败的证据。[Anthropic官方指标说明](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/data/docs/INSILICO.md)同样区分自洽性与准确性。

## 6. 能够确认并保留的内容

- PR #1列出的11项输入/参考库哈希全部匹配。
- 19条原序列独立IMGT编号末位均为127；追加S后的19条末位均为128。因此对本次h-NbBCII10来源构建，补S并重新评分有依据。不要把某个固定FR4正则推广成所有VHH框架的完整性定义。[IMGT区域定义](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html)。
- PR #2阈值表可以复算：公开数据固定seed0、ef2full、设计靶标与设计化学计量，剔除120条实验结论不确定的Mature GDF-8，剩1320条，其中354阳性、966阴性。阈值0.6得到TP309、FP535、FN45、TN431，召回87.29%、阳性预测值36.61%。
- 保持既定0.6、不为凑数量放宽门槛是合理执行选择；该跨靶点参考集没有TNFR1，不能据此宣称TNFR1的最优阈值或成功率已获验证。官方文档明确不跨靶点、化学计量及预测器直接比较ipSAE。
- 单种子和五种子最大值仍需分别报告。在有效1320条中，从seed0改为五种子取最大，有73条改变0.6过线状态（13个阳性、60个阴性）。不能仅凭中位差小将两者等同。本次seed0复算差值中位0.00922、均值0.03616；PR所写0.012/0.050没有附生成脚本、逐行选种与版本记录，未逐项复现。
- BC2在当前配置中的序列初筛率15/250高于RFD4的4/520及RFantibody的0/357；这支持后续资源比较，不支持直接预测每GPU小时收益或完全放弃探索。

公开数据来源与SHA-256记录在`results/calibration_recomputed.json`，复算脚本为`scripts/recompute.py`。原始文件为Anthropic发布的`data/tables/design_summary.parquet`和`data/tables/insilico/cofold_predictions.parquet`。[实验标签说明](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/data/docs/WETLAB.md)。

两份PR的GitHub CI通过，但其工作流只运行原有冻结结果验证，没有覆盖新增合规工具和外部校准统计；CI通过不能代替上述审查。

## 7. 当前扩量方案自身的缺口

只用h-NbBCII10是执行者为复用已验证原生输入、控制变量所作的实施选择，不是用户限制，也不是VHH模型必须单框架。该选择用于小批链路验证有明确用途；扩到约5000条时仍沿用一个框架，框架多样性没有同步扩大。

不仅生成输入只有一个框架，当前提取器还硬编码该框架的保守核心，其他框架输出会被拒收。多框架扩量需要同时调整原生结构条件、每个框架的元数据与编号校验、采样配额和最终结构聚类，不能只换一份FASTA。

建议后续修订：保留跨路线对照，加入经过完整性验证的多个VHH框架；序列数、独立亲本数及CDR/结合姿态结构簇数分开报告；按新增独立过线亲本/GPU小时评估采样配比。当前已提交扩量作业未因本次审查被停止或改配。


## 8. 上传内容和复现

这是已经完成的CPU独立复核结果；GPU扩量仍以其单独的PR #3和实时Slurm状态为准。本目录包含逐条IMGT/CDR结果、原始编号记录、逐设计seed0分数与实验标签、阈值表、路线/亲本计数、错误复现案例、输入和输出哈希。`sources/`是同伴PR的固定版本取证副本，含已发现问题的旧检查器，仅作为回归案例使用。

Windows原检查器产生的复现CSV使用系统默认编码；`results/`中提供UTF-8版本，原始字节保存在`provenance/original_novelty_backend_repro_cp936.csv`。

在仓库根目录：

```bash
python analyses/20261106-partner-pr-independent-review/scripts/validate_review.py
python analyses/20261106-partner-pr-independent-review/scripts/recompute.py --local-only
```

第一条仅用标准库，从原始编号和逐设计CSV复算关键结论并核对所有文件哈希；第二条重现旧检查器的CDR2误判和未运行新颖性后端却返回PASS的问题。文件写到忽略的`artifacts/`目录，原始证据保持不变。

若要从外部Parquet重新做完整统计，先安装pandas与pyarrow，再运行：

```bash
python analyses/20261106-partner-pr-independent-review/scripts/fetch_public_data.py
python analyses/20261106-partner-pr-independent-review/scripts/recompute.py
```

下载使用固定Hugging Face修订版本并核对SHA-256；数据和依赖不提交到Git。ANARCII原始Slurm执行脚本和资源设置保存在`provenance/executed/`，模型推理只能在Slurm计算节点执行。
