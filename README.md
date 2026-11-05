# Binder Design 协作项目

本仓库用于蛋白 binder 设计的代码、计算结果和团队协作。当前首个项目为 **TNFR1 单域 VHH 设计**。

**当前结果属于计算筛选候选，尚未完成结合实验或最终提交审计。**

| 路线 | 已统一评分 | 计算初筛通过 |
| --- | ---: | ---: |
| RFD4 | 520 | 4 |
| RFantibody | 357 | 0 |
| BindCraft2 | 250 | 15 |
| 合计 | **1,127** | **19** |

19 条初筛序列来自 8 条父设计轨迹。当前最终可提交数为 0；尚需完成公开序列原创性、完整靶标糖链可及性、多样性和最终选择。GPU 配额消耗为 74.20 小时，其中包含已记录的 BC2 故障开销。

## 从哪里开始

- [项目结果与目录说明](campaigns/tnfr1-vhh-20261103/README.md)
- [全部候选评分表](campaigns/tnfr1-vhh-20261103/results/all_candidates.csv)
- [19 条初筛序列](campaigns/tnfr1-vhh-20261103/shortlist/candidates.fasta)
- [19 条候选的结构](campaigns/tnfr1-vhh-20261103/shortlist/structures)
- [冻结评分协议](campaigns/tnfr1-vhh-20261103/config/scoring_protocol_v1.json)
- [工作流程与待完成工作](docs/WORKFLOW.md)
- [超算使用与复现边界](docs/HPC.md)
- [协作约定](CONTRIBUTING.md)

## 克隆和协作

仓库保持私有，协作者需要获得仓库访问权限后才能克隆或推送。

```bash
git clone https://github.com/RFjason/binder-design.git
cd binder-design
git switch -c codex/my-analysis
python tools/validate_results.py
# 修改分析代码或新增分析文件
git add <具体文件>
git commit -m "Add candidate analysis"
git push -u origin codex/my-analysis
```

通过 Pull Request 合并到 `main`。更新本地前先提交或保存自己的改动，再执行 `git pull --ff-only`。

## 原始计算附件

[TNFR1 初始结果 Release](https://github.com/RFjason/binder-design/releases/tag/tnfr1-v0.1.0) 保存所有 1,127 条候选的原始生成结构、复合物与单体预测、PAE 和输入映射。Git 内保留全部序列与评分，以及 19 条初筛候选的结构。

```bash
gh release download tnfr1-v0.1.0 --repo RFjason/binder-design \
  --pattern 'tnfr1-full-predictions-part*.zip' --dir artifacts/tnfr1-v0.1.0
python tools/validate_results.py --archives artifacts/tnfr1-v0.1.0
```

附件的 SHA-256 和目录说明位于项目的 `provenance/archive_manifest.json`。解压各分卷到同一个本地目录；它们是独立 ZIP 文件，不是需要拼接的分割压缩文件。

本仓库不包含模型权重、SSH 私钥、访问令牌或超算登录凭证。
