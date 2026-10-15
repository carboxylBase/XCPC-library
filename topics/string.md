# 字符串

## 变换不变量与非空可达类

- [QOJ18970 - abc - 变换](../training/2026/10/2026-10-07/QOJ18970-ABC-Transformation/notes.md)
  - a/b/c 映射 1/2/3 后展开与抵消都保持异或和；还需证明非空同类串之间的可达性，空串不能直接展开。
  - 字符编码不是四状态数组下标；合成第三种字符后需要减 'a'，或直接用 3-j-x。

## 模数性质与子串

- [QOJ 19008 - 神话子序列](../training/2026/09/2026-09-03/QOJ19008-Mythical-Subsequence/notes.md)
  - 十进制数模 `9` 等于数位和模 `9`；所有子串均非零等价于模 `9` 前缀和互不相同，因此合法串长度最多为 `8`。
  - 逐个追加字符时，新增子串恰好是所有新后缀，可以维护后缀余数集合而不保存全部子串。

## 反串配对与后缀 LCP 求和

- [QOJ18967 - 字符串的 Porder 之和](../training/2026/10/2026-10-14/QOJ18967-Porder-Sum/notes.md)
  - 逐层贡献转等长反串片段配对；去掉位置顺序后用 (C+P)/2 校正，P 是回文出现次数。
  - 原串与反串混合建 SA；height 链的路径最小值表达任意两后缀 LCP，降序并查集求跨串配对和。
  - Kasai 线性求 height，Manacher 按中心累加回文个数；倍增加 sort 当前为 O(n log^2 n)。

## 隐式 Trie 与字典序遍历

- [QOJ15323 - Yet another mailbox problem](../training/2026/10/2026-10-15/QOJ15323-Yet-Another-Mailbox-Problem/notes.md)
  - 先输出完整前缀，再按下一权值升序处理扩展；达到 k 条即停止，图上的路径不必显式存储字符串。
