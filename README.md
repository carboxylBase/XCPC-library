# XCPC Library

carboxylBase 的 XCPC 算法模板与个人训练归档。

## 目录

- `cpp/`：原有的算法模板、比赛工作区与代码速查资料。
- `training/`：按训练日期保存题目、代码与复盘，并维护尚待补写实现的题目清单。
- `topics/`：按算法专题建立题目索引，不重复存放源代码。
- `template/`：经过整理和验证的比赛模板。
- `skills/`：供后续 agent 使用的仓库维护规则。

## 归档约定

题目目录使用 `平台 + 题号 + 英文短名`，例如：

```text
training/2026/08/2026-08-26/CF1559D2-Mocha-and-Diana/
```

每道题至少包含：

- `notes.md`：记录真实思考过程、关键观察、证明、错误和重做情况。

用户提供实现后，再归档：

- `solution.cpp`：清理过调试代码、可以独立编译提交的最终实现。

已经口头推出完整思路、但尚未提供代码的题目，统一记录在 [`training/PENDING.md`](training/PENDING.md)；代码归档并验证后从清单移除。

维护训练记录时遵循 [`skills/training-log/SKILL.md`](skills/training-log/SKILL.md)。
