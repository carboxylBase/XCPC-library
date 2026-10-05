# 图论

## 并查集

- [CF1559D2 - Mocha and Diana](../training/2026/08/2026-08-26/CF1559D2-Mocha-and-Diana/notes.md)
  - 两个森林同时加边；用连通块数量建立答案上界。
  - 固定锚点 `1`，先合并共同的外部块，再将互补类型的连通块代表配对。

- [CF506D - Mr. Kitayuta's Colorful Graph](../training/2026/08/2026-08-31/CF506D-Mr-Kitayutas-Colorful-Graph/notes.md)
  - 按颜色建立单色连通块；轻颜色枚举连通块内点对，重颜色逐个检查询问。
  - 只记录实际询问点对，将轻颜色贡献表的空间控制在 `O(q)`。

## 模乘排列与函数约束

- [CF603B - Moodular Arithmetic](../training/2026/10/2026-10-05/CF603B-Moodular-Arithmetic/notes.md)
  - 固定非零 k 时，x -> kx mod p 是排列；并查集合并相邻点后，每个连通块就是一个环。
  - 起点函数值有 p 种选择，沿环传播自动闭合；k=0 的不可逆约束单独处理。
