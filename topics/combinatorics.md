# 组合计数

## 格路与组合数

- [CF559C - Gerald and Giant Chess](../training/2026/09/2026-09-01/CF559C-Gerald-and-Giant-Chess/notes.md)
  - 单调格路数用组合数计算；预处理阶乘与逆阶乘后单次查询为 `O(1)`。
  - 将终点作为额外普通点，与障碍点共同参与计数 DP。

## 小下标、大上标组合数

- [CF451E - Devu and Flowers](../training/2026/09/2026-09-02/CF451E-Devu-and-Flowers/notes.md)
  - 当组合数上标可达 `1e14`、下标却至多为 19 时，直接计算短下降幂并乘小阶乘的逆元。
  - 广义组合数满足 `C(-n,t)=(-1)^t C(n+t-1,t)`，可用于解释负整数次幂的展开。
