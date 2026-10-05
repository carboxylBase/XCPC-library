# 数论

## 阶、原根与解数统计

- [QOJ20240 - Exponent](../training/2026/10/2026-10-05/QOJ20240-Exponent/notes.md)
  - 欧拉定理配合带余除法得出阶整除 phi(n)；有原根时 a^d=1 的解数为 gcd(d,phi(n))。
  - CRT 使各质数幂的解数相乘；2^t 的奇数表示为 (-1)^b 5^k，t>=3 时计数为 gcd(d,2)*gcd(d,2^(t-2))。
  - 全局约数与局部 phi 不能混淆；约数莫比乌斯变换恢复精确阶分布。

## 模乘排列与函数约束

- [CF603B - Moodular Arithmetic](../training/2026/10/2026-10-05/CF603B-Moodular-Arithmetic/notes.md)
  - 固定 k 下乘 k 构成非零余数排列，每圈起点函数值有 p 种选择；可以直接并查集找环，不必求原根。
  - 区分固定 k 与对所有 k 成立；单独处理 k=0 和 k=1。

## 离散对数与线性同余

- [CF1106F - Lunar New Year and a Recursive Sequence](../training/2026/10/2026-10-05/CF1106F-Lunar-New-Year-and-a-Recursive-Sequence/notes.md)
  - 原根将乘法递推转成指数线性递推，指数模 p-1。
  - 矩阵求系数 A，BSGS 将输入 m 转成指数 z，exgcd 解 A*x=z (mod p-1)；三段任务不能互相替代。
