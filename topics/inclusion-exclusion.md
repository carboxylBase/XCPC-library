# 容斥

## 按唯一首个坏事件分类

- [CF559C - Gerald and Giant Chess](../training/2026/09/2026-09-01/CF559C-Gerald-and-Giant-Chess/notes.md)
  - 不显式枚举黑格子集，而是将每条非法路径归入它遇到的第一个黑格。
  - 排序特殊点后做 `O(n^2)` 计数 DP，避免普通容斥中的重复扣除。

## 有上界的非负整数解

- [CF451E - Devu and Flowers](../training/2026/09/2026-09-02/CF451E-Devu-and-Flowers/notes.md)
  - 将每个上界写进因子 `1-x^{f_i+1}`，枚举子集后，符号 `(-1)^|T|` 就是对超出上界事件的容斥。
  - `n <= 20` 时直接展开分子的 `2^n` 项，剩余次数用插板法计数。
