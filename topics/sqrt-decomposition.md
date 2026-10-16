# 根号分治

## 按集合大小分治

- [CF506D - Mr. Kitayuta's Colorful Graph](../training/2026/08/2026-08-31/CF506D-Mr-Kitayutas-Colorful-Graph/notes.md)
  - 一个含 `k` 条边的轻颜色用 `O(k^2)` 枚举连通点对，利用 `sum(k^2) <= Bm` 控制总代价。
  - 重颜色至多 `m/B` 种，逐个扫描询问，平衡得到 `O(mB + qm/B)`。
  - 统一询问点对方向并使用哈希表可以进一步降低常数；本次使用有序映射的实现也在时限内通过。

## 按值域二进制位折半

- [QOJ15322 - Round screws](../training/2026/10/2026-10-16/QOJ15322-Round-Screws/notes.md)
  - 18 位拆为高低各 9 位；按历史高位分组，表内预存所有查询低位的最小值，更新固定一行、查询固定一列，各 512 次。

## 分块维护前缀极值

- [QOJ15324 - No more regrets](../training/2026/10/2026-10-16/QOJ15324-No-More-Regrets/notes.md)
  - 整块维护前缀最小值、最大值及其乘积的前缀和；区间加和赋值用懒标记，散块修改后重建。
  - 查询通过两个单调分界分成三段；用户最终 B=1500 版本 AC，注意重建次数与常数。
