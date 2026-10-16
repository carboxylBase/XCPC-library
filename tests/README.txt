Template regression tests (2026-10-16)

Run: python tests/template_regression.py
Requires Python 3 and GCC with C++17 and __int128. Set CXX to the compiler path if needed.
Tests use temporary copies, substituting all.hpp with bits/stdc++.h; repository files are not modified.
Scope: syntax-check template .cpp files excluding WorkSpace, plus 15 targeted regression/random test groups.
Passing syntax checks is not proof of algorithm correctness for untested templates.

Alternative names after separating duplicate declarations:
Lca / LcaBinary (Euler RMQ / binary lifting)
FHQ_Treap / FHQ_TreapValue
KD_Tree / KD_TreeDistance
LCT / LCTSubtree / LCTAffine
LiChao / LiChaoFloat
ChthollyTree / IntervalSet
DynamicSegTree / DynamicSegTreePoint; merge/split namespace SegTreeMergeSplit
SAM / SAMLowercase; array trie namespace TrieBits, vector class Trie
The first variant's name is retained except where a namespace/type collision required separation.

Matrix modulus: Mat<long long> defaults to 1000000007; Mat<long long,998244353> selects another modulus.
LinearBasis::get_mn returns the minimum nonempty subset XOR, or -1 if no element has been inserted.
Marked framework templates still require problem-specific aggregates and lazy-tag operations.
