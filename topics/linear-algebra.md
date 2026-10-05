# 线性代数

## 矩阵快速幂与线性递推

- [CF1106F - Lunar New Year and a Recursive Sequence](../training/2026/10/2026-10-05/CF1106F-Lunar-New-Year-and-a-Recursive-Sequence/notes.md)
  - 乘法递推先用原根转成指数线性递推；初始指数只有最后一项未知。
  - 状态升序排列时取 T^(n-k) 的右下角作为未知项系数；递推系数 b_i 不能误作初始值。
  - 矩阵计算模 p-1，最后还原函数值才模 p。
