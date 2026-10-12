# 有序集合

## 批次前驱查询与反悔候选

- [Treasure Island](../training/2026/09/2026-09-07/Unknown-Treasure-Island/notes.md)
  - 一个 multiset 按代价枚举反悔候选，另一个按资源出现时间维护批次，通过 lower_bound 的前驱定位严格早于战斗的最新金币。
  - 按批次扣除大量相同资源，结合候选的永久删除，以摊还分析保证 O(n log n)。
