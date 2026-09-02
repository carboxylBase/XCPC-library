# 生成函数

## 有限几何级数与系数提取

- [CF451E - Devu and Flowers](../training/2026/09/2026-09-02/CF451E-Devu-and-Flowers/notes.md)
  - 有上界变量用 `1+x+...+x^f=(1-x^{f+1})/(1-x)` 表示。
  - 分子通过子集枚举展开，分母使用 `[x^t](1-x)^{-n}=C(t+n-1,n-1)` 提取系数。
