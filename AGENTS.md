# 协作与超算规则

- 不在 sylogin/login 节点运行计算工作负载。登录节点仅用于轻量文件操作、传输、Slurm 作业管理及日志检查。
- RFD4、RFdiffusion、ProteinMPNN、AF2、RF2、ESMFold2、PyTorch/JAX 推理和批量结构分析均必须通过 Slurm 在计算节点运行。
- 不提交凭证、私钥、令牌、`.env`、模型权重或个人 SSH 配置。
- `campaigns/tnfr1-vhh-20261103` 是计算快照。不要直接改写其原始结果以提高分数；新实验或重新评分应保存到新目录，保留亲本和配置来源。
- 本项目成熟 VHH 限制为100–150 aa、20种标准天然氨基酸、完整框架和CDR，不含信号肽、标签或Fc融合。
- 主评分配置冻结为ESMFold2 Full、seed0、1 sample、无模板、无初始猜测、仅靶标MSA。变更评分仪器或参数必须作为新实验明确记录。
- `computational_quality_pass` 表示计算初筛。不能写成实验结合成功、已验证亲和力或最终可提交。
- 当前原创性、完整靶标糖链可及性和最终多样性审计未完成。未经完成的检查保留未完成状态。
- 每项新Slurm作业时限必须小于48小时，优先拆成短批次。启动设计前强制确认所需GPU后端可用，禁止静默回退CPU。
- 分支默认使用`codex/`前缀。合并前运行`python tools/validate_results.py`。
