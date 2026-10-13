# 有序集合

## 批次前驱查询与反悔候选

- [Treasure Island](../training/2026/09/2026-09-07/Unknown-Treasure-Island/notes.md)
  - 一个 multiset 按代价枚举反悔候选，另一个按资源出现时间维护批次，通过 lower_bound 的前驱定位严格早于战斗的最新金币。
  - 按批次扣除大量相同资源，结合候选的永久删除，以摊还分析保证 O(n log n)。

## multiset 维护圆方树方点权值

- [CF487E - Tourists](../training/2026/10/2026-10-12/CF487E-Tourists/notes.md)
  - 一个方点只存儿子圆点权值，begin() 为最小值；修改用 erase(find(old)) 删除一个实例，再 insert(new)，根不访问父方点集合。

- [QOJ18950 - 建设高铁](../training/2026/10/2026-10-13/QOJ18950-Railway-Construction/notes.md)
  - Kruskal 重构树将连通阈值变成集合 LCA；只需选中城市 DFS 序的最小、最大点。
  - 离线收集繁荣度/编号键，线段树维护有效数量与 DFS 序极值，按数量下降求前 k 名。用户确认 AC。

## 环形滑动窗口的最大值

- [CF106695L - 贪吃蛇](../training/2026/10/2026-10-13/CF106695L-Snake/notes.md)
  - multiset 保存窗口各方向的最远距离平方；rbegin() 取最大值，erase(find(value)) 只删除一个重复值。
