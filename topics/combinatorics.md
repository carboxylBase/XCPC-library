# 组合计数

## 格路与组合数

- [CF559C - Gerald and Giant Chess](../training/2026/09/2026-09-01/CF559C-Gerald-and-Giant-Chess/notes.md)
  - 单调格路数用组合数计算；预处理阶乘与逆阶乘后单次查询为 `O(1)`。
  - 将终点作为额外普通点，与障碍点共同参与计数 DP。

## 小下标、大上标组合数

- [CF451E - Devu and Flowers](../training/2026/09/2026-09-02/CF451E-Devu-and-Flowers/notes.md)
  - 当组合数上标可达 `1e14`、下标却至多为 19 时，直接计算短下降幂并乘小阶乘的逆元。
  - 广义组合数满足 `C(-n,t)=(-1)^t C(n+t-1,t)`，可用于解释负整数次幂的展开。

## 满射计数

- [CF1342E - Placing Rooks](../training/2026/09/2026-09-02/CF1342E-Placing-Rooks/notes.md)
  - 从 `n` 个有标号元素到 `m` 个有标号集合的满射数为 `sum(-1)^(m-i) C(m,i)i^n`。
  - 先由结构确定非空集合数 `m=n-k`，再选出具体集合并计算满射。

## 独立分量的函数计数

- [CF603B - Moodular Arithmetic](../training/2026/10/2026-10-05/CF603B-Moodular-Arithmetic/notes.md)
  - 函数不要求单射，每个模乘环可自由选择起点函数值，包括 0；C 个非零环贡献 p^C。
- [QOJ20240 - Exponent](../training/2026/10/2026-10-05/QOJ20240-Exponent/notes.md)
  - CRT 让局部方程的解数相乘；先统计阶整除 d 的数量，再反演得到阶恰好为 d 的数量并加权。
