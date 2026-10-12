#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

/*
1 ~ n 是圆点
n + 1 ~ tot 是方点
*/

struct BlockCutTree {
    int n, tot, timer = 0;
    vector<vector<pair<int, int>>> g;
    vector<vector<int>> tree;
    vector<int> dfn, low, stk;
    int edgeCnt = 0;

    explicit BlockCutTree(int n)
        : n(n), tot(n), g(n + 1), tree(2 * n + 1),
          dfn(n + 1), low(n + 1) {}

    void addEdge(int u, int v) {
        if (u == v) return;
        int id = edgeCnt++;
        g[u].push_back({v, id});
        g[v].push_back({u, id});
    }

    void link(int u, int v) {
        tree[u].push_back(v);
        tree[v].push_back(u);
    }

    void dfs(int u, int parentEdge) {
        dfn[u] = low[u] = ++timer;
        stk.push_back(u);

        for (auto [v, id] : g[u]) {
            if (id == parentEdge) continue;

            if (!dfn[v]) {
                dfs(v, id);
                low[u] = min(low[u], low[v]);

                if (low[v] >= dfn[u]) {
                    int b = ++tot;
                    link(b, u);

                    while (true) {
                        int x = stk.back();
                        stk.pop_back();
                        link(b, x);
                        if (x == v) break;
                    }
                }
            } else {
                low[u] = min(low[u], dfn[v]);
            }
        }
    }

    void build() {
        for (int u = 1; u <= n; ++u) {
            if (!dfn[u]) {
                dfs(u, -1);
                stk.pop_back();
            }
        }
        tree.resize(tot + 1);
    }
};

#define LS (rt << 1)
#define RS (rt << 1 | 1)

struct SegTree {
    struct Node {
        int mn;
        Node() {
            mn = INT_MAX;
        }
    } nodes[N * 4];

    Node merge(Node L, Node R) {
        Node M;
        M.mn = min(L.mn, R.mn);
        return M;
    }

    void build(int rt, int l, int r) {
        if (l == r) {
            return;
        }
        int mid = (l + r) >> 1;
        build(LS, l, mid), build(RS, mid + 1, r);
        nodes[rt] = merge(nodes[LS], nodes[RS]);
    }

    void pd(int rt) {

    }

    void update(int rt, int l, int r, int ql, int qr, int val) {
        if (ql > qr || l > qr || ql > r) return;
        if (ql <= l && r <= qr) {
            return;
        }
        int mid = (l + r) >> 1;
        pd(rt);
        if (ql <= mid) {
            update(LS, l, mid, ql, qr, val);
        }
        if (qr >= mid + 1) {
            update(RS, mid + 1, r, ql, qr, val);
        }
        nodes[rt] = merge(nodes[LS], nodes[RS]);
        return;
    }

    void modify(int rt, int l, int r, int q, int val) {
        if (l > q || r < q) return;
        if (l == r) {
            nodes[rt].mn = val;
            return;
        }
        int mid = (l + r) >> 1;
        pd(rt);
        if (q <= mid) {
            modify(LS, l, mid, q, val);
        }
        if (q >= mid + 1) {
            modify(RS, mid + 1, r, q, val);
        }
        nodes[rt] = merge(nodes[LS], nodes[RS]);
        return;
    }

    Node query(int rt, int l, int r, int ql, int qr) {
        if (ql > qr) {
            return Node();
        }
        if (ql <= l && r <= qr) {
            return nodes[rt];
        }
        int mid = (l + r) >> 1;
        pd(rt);
        if (ql > mid) {
            return query(RS, mid + 1, r, ql, qr);
        } else if (qr < mid + 1) {
            return query(LS, l, mid, ql, qr);
        } else {
            return merge(
                query(LS, l, mid, ql, qr),
                query(RS, mid + 1, r, ql, qr)
            );
        }
    }
};

struct TreeChainSeg {
    vector<vector<int>> g;
    SegTree seg;
    vector<int> fa, top, siz, son, dep, dfn, id;
    int cnt, n, cnt2;

    void dfs1(int u) {
        siz[u] = 1;
        for (auto v : g[u]) {
            if (fa[u] == v) continue;
            dep[v] = dep[u] + 1;
            fa[v] = u;
            dfs1(v);
            siz[u] += siz[v];
            if (!son[u] || siz[son[u]] < siz[v]) {
                son[u] = v;
            }
        }
        return;
    }

    void dfs2(int u, int tp) {
        dfn[u] = ++cnt;
        id[cnt] = u;
        top[u] = tp;
        if (son[u]) {
            dfs2(son[u], tp);
        }
        for (auto v : g[u]) {
            if (v == fa[u] || v == son[u]) continue;
            dfs2(v, v);
        }
        return;
    }

    void init(vector<vector<int>>& G) {
        g = G;
        n = (int)g.size() - 1;
        cnt = 0, cnt2 = n + 1;
        fa.clear(); fa.resize(n + 1, 0);
        top.clear(); top.resize(n + 1, 0);
        siz.clear(); siz.resize(n + 1, 0);
        son.clear(); son.resize(n + 1, 0);
        dep.clear(); dep.resize(n + 1, 0);
        dfn.clear(); dfn.resize(n + 1, 0);
        id.clear(); id.resize(n + 1, 0);

        dfs1(1), dfs2(1, 1);
        seg.build(1, 1, n);
        return;
    }

    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                u = fa[top[u]];
            } else {
                v = fa[top[v]];
            }
        }
        if (dep[u] > dep[v]) {
            return v;
        }
        return u;
    }

    pii query(int u, int v) {
        int ans = INT_MAX;
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                ans = min(
                    ans,
                    seg.query(1, 1, n, dfn[top[u]], dfn[u]).mn
                );
                u = fa[top[u]];
            } else {
                ans = min(
                    ans,
                    seg.query(1, 1, n, dfn[top[v]], dfn[v]).mn
                );
                v = fa[top[v]];
            }
        }
        if (dep[u] > dep[v]) {
            swap(u, v);
        }
        ans = min(ans, seg.query(1, 1, n, dfn[u], dfn[v]).mn);
        return pii(ans, u);
    }
} seg;

void solve() {
    int n, m, q;
    cin >> n >> m >> q;

    vector<int> W(n + 1, 0);
    for (int i = 1; i < n + 1; i++) cin >> W[i];

    BlockCutTree solver(n);
    for (int i = 1; i < m + 1; i++) {
        int u, v;
        cin >> u >> v;
        solver.addEdge(u, v);
    }
    solver.build();

    seg.init(solver.tree);
    vector<int> fa(solver.tot + 1, 0);

    auto dfs = [&](auto&& self, int u) -> void {
        for (auto v : solver.tree[u]) {
            if (v == fa[u]) {
                continue;
            }
            fa[v] = u;
            self(self, v);
        }
        return;
    };
    dfs(dfs, 1);

    vector<multiset<int>> vals(solver.tot - n);
    for (int i = 1; i <= n; i++) {
        if (fa[i] > n) {
            vals[fa[i] - n - 1].insert(W[i]);
        }
    }

    for (int i = 1; i <= solver.tot; i++) {
        if (i <= n) {
            seg.seg.modify(1, 1, solver.tot, seg.dfn[i], W[i]);
        } else {
            seg.seg.modify(
                1, 1, solver.tot, seg.dfn[i],
                *vals[i - n - 1].begin()
            );
        }
    }

    while (q--) {
        char opt;
        cin >> opt;

        if (opt == 'C') {
            ll a, w;
            cin >> a >> w;

            seg.seg.modify(1, 1, solver.tot, seg.dfn[a], w);

            if (fa[a] > n) {
                vals[fa[a] - n - 1].erase(
                    vals[fa[a] - n - 1].find(W[a])
                );
                vals[fa[a] - n - 1].insert(w);

                seg.seg.modify(
                    1, 1, solver.tot, seg.dfn[fa[a]],
                    *vals[fa[a] - n - 1].begin()
                );
            }

            W[a] = w;
        } else {
            ll a, b;
            cin >> a >> b;

            if (a == b) {
                cout << W[a] << '\n';
                continue;
            }

            pii res = seg.query(a, b);
            if (res.second > n) {
                res.first = min(res.first, W[fa[res.second]]);
            }
            cout << res.first << '\n';
        }
    }
    return;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
