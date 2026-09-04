# 构造

## 奇偶性与缩放

- [CF1270E - Divide Points](../training/2026/09/2026-09-03/CF1270E-Divide-Points/notes.md)
  - 先用平方距离模 2、模 4 区分同组与异组点对；所有点同奇偶时统一除以 2，逐层提取坐标差的最低有效二进制位。

## 倍增与覆盖

- [CF1408F - Two Different](../training/2026/09/2026-09-05/CF1408F-Two-Different/notes.md)
  - 让两个内部相等的等长块按同一方向逐位调用任意函数，可以倍增出更大的相等块。
  - 用前后两个最大 2 的幂长度的重叠区间覆盖任意长度，只留下至多两个最终值。

## 通信与编码

- [QOJ 18436 - Nailoong vs. Bombloong 2](../training/2026/08/2026-08-30/QOJ18436-Nailoong-vs-Bombloong-2/notes.md)
  - 将染色后的割边数差转化为节点儿子数，再按 BFS 分组重建同构树。
  - 利用固定信息和全局总和消去三个无需传递的量，是满足严格通信长度的通用思路。
