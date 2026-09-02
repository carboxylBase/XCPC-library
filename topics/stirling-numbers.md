# 斯特林数

## 第二类斯特林数与满射

- [CF1342E - Placing Rooks](../training/2026/09/2026-09-02/CF1342E-Placing-Rooks/notes.md)
  - `S(n,m)` 把 `n` 个有标号元素划分为 `m` 个非空无标号集合；给集合加上标号后，满射数为 `m!S(n,m)`。
  - 本题从指数生成函数与容斥独立推出 `m!S(n,m)=sum(-1)^(m-i)C(m,i)i^n`。
