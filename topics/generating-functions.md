# 生成函数

## 有限几何级数与系数提取

- [CF451E - Devu and Flowers](../training/2026/09/2026-09-02/CF451E-Devu-and-Flowers/notes.md)
  - 有上界变量用 `1+x+...+x^f=(1-x^{f+1})/(1-x)` 表示。
  - 分子通过子集枚举展开，分母使用 `[x^t](1-x)^{-n}=C(t+n-1,n-1)` 提取系数。

## 指数生成函数与满射

- [CF1342E - Placing Rooks](../training/2026/09/2026-09-02/CF1342E-Placing-Rooks/notes.md)
  - 非空有标号集合对应 `e^x-1`，把 `n` 个有标号元素分到 `m` 个非空有标号集合对应 `n![x^n](e^x-1)^m`。
  - 展开为 `e^{ix}` 后，可直接使用 `[x^n]e^{ix}=i^n/n!`，无需多项式乘法。
