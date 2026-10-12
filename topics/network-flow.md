# 网络流

## 二分图匹配与最小割

- [CF2026E - Best Subsequence](../training/2026/08/2026-08-31/CF2026E-Best-Subsequence/notes.md)
  - 将数字节点放在左部、二进制位节点放在右部，源点和汇点侧容量均为 `1`，最大流即最大匹配。
  - 冲突图的最大独立集通过 König 定理转化为最大匹配，最终答案为 `n - maxflow`。

