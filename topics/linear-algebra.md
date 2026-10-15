# 线性代数

## 矩阵快速幂与线性递推

- [CF1106F - Lunar New Year and a Recursive Sequence](../training/2026/10/2026-10-05/CF1106F-Lunar-New-Year-and-a-Recursive-Sequence/notes.md)
  - 乘法递推先用原根转成指数线性递推；初始指数只有最后一项未知。
  - 状态升序排列时取 T^(n-k) 的右下角作为未知项系数；递推系数 b_i 不能误作初始值。
  - 矩阵计算模 p-1，最后还原函数值才模 p。

## 线性基与掩码投影

- [QOJ15320 - Gemcrate](../training/2026/10/2026-10-15/QOJ15320-Gemcrate/notes.md)
  - 固定掩码的按位 AND 是异或空间的线性映射；原基投影后重新建基，检查能否覆盖要求的 1 位。
  - 最终 AND 为 0 的位置不等于单组 XOR 必须为 0；[11,13] 是反例。临时变量必须保持 64 位。
