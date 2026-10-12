# 记号约定

+ ACAM：AC 自动机。
+ Trie：字典树。
+ $S$：文本串。

# fail 树

考虑如何在文本串中找出所有能够与模式串匹配的子串终点。首先，将所有模式串插入 Trie，并标记所有作为某个模式串终点的节点。

接下来考虑匹配过程。假设当前位于节点 $u$，准备匹配文本串 $S$ 的第 $i$ 个字符 $c$。如果 $u$ 存在字符 $c$ 对应的儿子，那么直接转移到这个儿子即可；否则，需要找到 $u$ 的最长真后缀，使得这个后缀存在字符 $c$ 对应的儿子。于是，$u$ 需要不断跳向自己的真后缀。由节点指向其最长真后缀的边称为 fail 边，所有 fail 边构成 fail 树。

考虑 fail 树的构建方式。先建好 Trie，然后自顶向下更新 fail 数组。假设当前位于节点 $u$，考虑 $u$ 经过字符 $c$ 到达的儿子 $v$：

+ 若 $v$ 存在，则

$$
\mathrm{fail}[v]=\mathrm{ch}[\mathrm{fail}[u]][c].
$$

+ 若 $v$ 不存在，则

$$
\mathrm{ch}[u][c]=\mathrm{ch}[\mathrm{fail}[u]][c].
$$

需要注意的是，此时 $\mathrm{ch}$ 的语义已经发生变化：它不再只表示 Trie 中从 $u$ 经过字符 $c$ 到达的儿子，而是表示在 AC 自动机上，从当前状态 $u$ 接受字符 $c$ 后能够到达的状态。换句话说，当 Trie 中不存在对应儿子时，$\mathrm{ch}[u][c]$ 会被补成沿 fail 边回退后得到的转移。

这样一来，在文本串匹配过程中，每次接受一个新的字符，都可以从上一个匹配到的最长后缀，转移到新的最长匹配后缀。然而，仅靠这个转移还不足以得到所有模式串终点对应的位置。

例如，存在模式串 `abc` 和 `bc`，文本串为 `abcd`。当当前位置已经匹配到 `abc` 时，节点 `abc` 本身对应一个模式串终点；同时，`bc` 也是 `abc` 的后缀，也应该在当前位置被统计到。如果只记录当前最长匹配后缀，就会漏掉这些较短的模式串。因此，需要从 `abc` 对应的节点沿 fail 树向上跳，找出所有作为 `abc` 后缀的模式串节点，然后再回到 `abc` 对应的节点，继续接受后面的字符 `d`。

但是，暴力跳 fail 树的复杂度显然不可控。为此，额外维护 $\mathrm{out}[u]$，表示节点 $u$ 的最长模式串后缀。转移如下：

+ 若 $\mathrm{fail}[u]$ 是模式串终点，则

$$
\mathrm{out}[u]=\mathrm{fail}[u].
$$

+ 若 $\mathrm{fail}[u]$ 不是模式串终点，则

$$
\mathrm{out}[u]=\mathrm{out}[\mathrm{fail}[u]].
$$

上面两个递推式都要求自顶向下更新。参考代码如下：

```cpp
namespace ACAM {
    struct Node {
        int ch[26], idx, fail, out;
        Node() {
            for (int i = 0;i<26;i++) ch[i] = -1;
            idx = fail = out = -1;
        }
    } nodes[N];

    int rt, cnt_node;

    int newNode() {
        nodes[cnt_node] = Node();
        return cnt_node ++;
    }

    void init() {
        cnt_node = 0;
        rt = newNode();
        return;
    }

    // 这里默认使用 ACAM::rt, dep = 0
    void insert(int rt, int dep, string& s, int idx) {
        if (s.size() == dep) {
            nodes[rt].idx = idx;
            return;
        }
        if (nodes[rt].ch[s[dep] - 'a'] == -1) nodes[rt].ch[s[dep] - 'a'] = newNode();
        insert(nodes[rt].ch[s[dep] - 'a'], dep + 1, s, idx);
        return;
    }

    void process() {
        queue<int> q;
        
        for (int i = 0;i<26;i++) {
            int v = nodes[rt].ch[i];
            if (v != -1) {
                q.push(v);
                nodes[v].fail = rt;
            } else {
                nodes[rt].ch[i] = rt;
            }
        }

        while (!q.empty()) {
            int u = q.front(); 
            q.pop();

            for (int i = 0;i<26;i++) {
                int v = nodes[u].ch[i];
                if (v == -1) {
                    nodes[u].ch[i] = nodes[nodes[u].fail].ch[i];
                } else {
                    q.push(v);
                    nodes[v].fail = nodes[nodes[u].fail].ch[i];
                    if (nodes[nodes[v].fail].idx != -1) {
                        nodes[v].out = nodes[v].fail;
                    } else {
                        nodes[v].out = nodes[nodes[v].fail].out;
                    }
                }
            }
        }

        return;
    }
};
```

# 习题推荐

[$2026$ 牛客暑期多校 $L$ 题](https://ac.nowcoder.com/acm/contest/133876/L)

本题核心在于找出所有模式串在文本串中出现的问题