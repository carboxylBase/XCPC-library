# 前缀和

## 二维区域中的合法点权和

- [CF2275H - A Problem to Warm Up the Eyebrows](../training/2026/10/2026-10-23/CF2275H-A-Problem-to-Warm-Up-the-Eyebrows/notes.md)
  - 先按坐标奇偶筛选点权，再用左上前缀和与右上前缀和统计两种对角布局，将点对求和降为 O(nm)。
  - 复用数组前重置每格底值；第二遍查询 b[i-1][j+1]，同时排除同行与同列的重复贡献。
