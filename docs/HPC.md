# 超算与复现说明

## 安全规则

登录节点仅进行轻量文件操作、传输、Slurm管理和日志检查。所有设计、预测、批量结构分析和压缩导出均通过Slurm计算节点执行。

每个作业时限小于48小时；已有批次采用20–45分钟评分或30分钟生成。GPU后端初始化失败必须立即停止，不能在已申请GPU的作业里静默回退CPU。

## 历史脚本快照

`campaigns/tnfr1-vhh-20261103/scripts/hpc_snapshot/`是实际执行脚本的快照，包含原HPC绝对路径、环境路径和历史作业号。它们用于追溯，不能在其他集群不经适配就直接运行。

尤其不要直接执行历史`submit_*`、`shorten_*`、recovery作业管理脚本：其中可能引用本轮已经结束的Slurm数组。新轮次应新建配置和提交脚本，使用`sbatch`返回的新作业ID。

原项目存储位置记录在冻结协议及`provenance/job_ledger.json`中。原始模型权重与已安装工具不随Git仓库分发；同伴需要自行准备被授权的环境和上游权重。

## 本轮故障记录

初次BC2批量启动时cuSPARSE加载失败，JAX回退CPU，24个作业达到2小时时限而没有设计输出。修复为优先加载BC2虚拟环境自身NVIDIA库，并设置`JAX_PLATFORMS=cuda`和显式GPU检查。修复后的40批完整完成、无运行失败或超时。

74.20 GPU小时包含上述故障占用的GPU配额，不是有效GPU推理时间。Slurm记账优先于只覆盖正常退出作业的脚本耗时。

## 原始附件目录

Release分卷解压后，每个候选具有：

```text
designs/<candidate_id>/original.cif 或 original.pdb
designs/<candidate_id>/designed_complex.pdb
predictions/<candidate_id>/complex/model.cif
predictions/<candidate_id>/complex/pae.npz
predictions/<candidate_id>/complex/summary.json
predictions/<candidate_id>/complex/input_and_token_mapping.json
predictions/<candidate_id>/monomer/...
```

`pae.npz`同时保存原始PAE、按CIF顺序重排的PAE和链/token映射。所有比较需要沿用这些映射，不能直接假定矩阵行号等于PDB残基编号。

复合物预测包含任务指定的1–211目标序列；自一致性DockQ只比较成熟胞外域与VHH。RFantibody部分原始设计只有主链，源结构缺少SG标为不可评估；二硫键最终在完整ESMFold2单体中检查。
