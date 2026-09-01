# 容斥

## 按唯一首个坏事件分类

- [CF559C - Gerald and Giant Chess](../training/2026/09/2026-09-01/CF559C-Gerald-and-Giant-Chess/notes.md)
  - 不显式枚举黑格子集，而是将每条非法路径归入它遇到的第一个黑格。
  - 排序特殊点后做 `O(n^2)` 计数 DP，避免普通容斥中的重复扣除。
