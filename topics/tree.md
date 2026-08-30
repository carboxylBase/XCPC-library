# 树

## 序列编码与重构

- [QOJ 18436 - Nailoong vs. Bombloong 2](../training/2026/08/2026-08-30/QOJ18436-Nailoong-vs-Bombloong-2/notes.md)
  - 选择叶子作为根，用 BFS 顺序中的儿子数序列唯一编码一棵有根树。
  - 根和末尾叶子的儿子数固定，再用儿子数总和补出倒数第二个节点，将通信次数压到 `n - 3`。
