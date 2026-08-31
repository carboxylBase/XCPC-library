# 根号分治

## 按集合大小分治

- [CF506D - Mr. Kitayuta's Colorful Graph](../training/2026/08/2026-08-31/CF506D-Mr-Kitayutas-Colorful-Graph/notes.md)
  - 一个含 `k` 条边的轻颜色用 `O(k^2)` 枚举连通点对，利用 `sum(k^2) <= Bm` 控制总代价。
  - 重颜色至多 `m/B` 种，逐个扫描询问，平衡得到 `O(mB + qm/B)`。
  - 实现时应规范化并去重询问，重颜色只扫描不同点对；轻颜色用哈希查询，避免额外的 `log q`。
