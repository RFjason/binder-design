# 当前工作流

任务与构建体核查 → 表位/框架/CDR约束 → 三路线生成 → 独立ESMFold2重预测 → 计算筛选 → 终审与优化。

## 科学输入

- 靶标：人TNFR1，UniProt P19438-1。
- 设计构建体：成熟胞外域30–211；主评分构建体按任务要求为1–211。
- 核对9S3I与UniProt编号，并用完整成熟胞外域建立设计背景。
- 两组候选热点：CRD2 L96/H98/L100；CRD3 H134/L140/Q142。它们是结构候选热点，不是经突变实验确认的能量热点。
- VHH使用来源明确的h-NbBCII10通用框架，固定83 aa框架核心，设计三个CDR片段及结合姿态。不是从零重新设计免疫球蛋白框架。

## 生成与评分

RFD4联合生成结构和序列；RFantibody串联RFdiffusion、ProteinMPNN及RF2；BindCraft2梯度优化CDR并进行序列重设计。

所有候选统一使用ESMFold2 Full、seed=0、1 sample、无模板、无初始猜测、仅靶标MSA。主分数是双向ipSAE的较小值，PAE与距离cutoff均为10 Å。另做独立单体预测和设计—重预测的自一致性DockQ。

基础计算门槛：ipSAE≥0.6、单体平均pLDDT≥70、DockQ≥0.23、完整单体二硫键几何通过。DockQ参考界面不可评估的候选保存空值且不能通过门槛。

## 当前完成与未完成

已完成1127条候选评分和19条初筛候选的局部CDR置信度、表面性质、潜在糖基化基序检查。

尚未完成：

- TNFR1同靶标已知阳性/阴性对照分离校准；
- 完整UniRef90/公开binder库原创性审计；
- Anthropic提示词中要求的LCP局部序列组成约束；
- 完整靶标糖链可及性审计；
- 结构聚类、多样性和最终选择；
- 围绕优选亲本的系统多轮定向优化；
- 表达、结合和功能实验验证。

当前19条是计算初筛候选，最终可提交数暂记0。

## 与Anthropic流程的关系

借鉴自主表位选择、多工具编排、独立重预测、来源追踪和预算控制。按用户任务改成VHH和单seed主评分，没有完整复现其七生成路线、多预测模型终排、常驻多代理监控、对照校准、完整原创性筛选及实验验证。

参考：[公开单靶标提示词](https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/blob/main/prompts/prompts/single_target/EGFR.md)、[技术报告](https://www-cdn.anthropic.com/30bf50e22a01388bb29bf077ee3f244531594b7a.pdf)。
