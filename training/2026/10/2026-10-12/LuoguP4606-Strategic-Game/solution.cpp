#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

/*
1 ~ n 是圆点
n + 1 ~ tot 是方点
你要查询 u, v 直接的必经之点, 本质上就是找 u 和 v 之间的圆点, 包括 u v 在内直接查找
*/

struct BlockCutTree {
    int n, tot, timer = 0;
    vector<vector<pair<int, int>>> g;  // {邻点, 无向边编号}
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
                    int b = ++tot;  // 新建方点
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

    // 加完所有边后调用一次
    void build() {
        for (int u = 1; u <= n; ++u) {
            if (!dfn[u]) {
                dfs(u, -1);
                stk.pop_back();  // 该连通分量的 DFS 根还在栈内
            }
        }
        tree.resize(tot + 1);
    }
};

void solve() {
    int n, m; cin >> n >> m;
    BlockCutTree solver(n);
    for (int i = 1; i < m + 1; i++) {
        int u, v; cin >> u >> v;
        solver.addEdge(u, v);
    }
    solver.build();

    int M = 22;
    vector<vector<int>> fa(M, vector<int>(solver.tot + 1, 0));
    vector<int> dep(solver.tot + 1, 0);
    vector<vector<pii>> g(solver.tot + 1);
    vector<int> dfn(solver.tot + 1, 0);
    int cnt_dfn = 0;

    auto dfs = [&](auto&& self, int u) -> void {
        dfn[u] = ++ cnt_dfn;
        for (auto v : solver.tree[u]) {
            if (v == fa[0][u]) {
                continue;
            }
            fa[0][v] = u;
            dep[v] = dep[u] + 1;
            self(self, v);
        }
        return;
    };
    dep[1] = 1;
    dfs(dfs, 1);

    for (int j = 1; j < M; j++) {
        for (int i = 1; i <= solver.tot; i++) {
            fa[j][i] = fa[j - 1][fa[j - 1][i]];
        }
    }

    vector<vector<int>> dp(M, vector<int>(solver.tot + 1, 0));
    for (int i = 1; i <= solver.tot; i++) {
        if (i <= n) {
            dp[0][i] = 1;
        }
    }
    for (int j = 1; j < M; j++) {
        for (int i = 1; i <= solver.tot; i++) {
            dp[j][i] = dp[j - 1][i] + dp[j - 1][fa[j - 1][i]];
        }
    }

    int q; cin >> q;
    while (q--) {
        int s; cin >> s;
        vector<int> b(s);
        for (auto &v : b) {
            cin >> v;
        }
        sort(b.begin(), b.end(), [&](int x, int y) {
            return dfn[x] < dfn[y];
        });
        auto lca = [&](int x, int y) -> int {
            for (int i = M - 1;i>=0;i--) {
                if (dep[x] > dep[y]) swap(x, y);
                if (dep[fa[i][y]] >= dep[x]) {
                    y = fa[i][y];
                }
            }
            if (x == y) return x;
            for (int i = M - 1;i>=0;i--) {
                if (fa[i][x] != fa[i][y]) {
                    x = fa[i][x], y = fa[i][y];
                }
            }
            return fa[0][x];
        };
        auto cal = [&](int x, int y) -> int {
            int res = 0;
            for (int i = M - 1;i>=0;i--) {
                if (dep[x] > dep[y]) swap(x, y);
                if (dep[fa[i][y]] >= dep[x]) {
                    res += dp[i][y];
                    y = fa[i][y];
                }
            }
            if (x == y) return res;
            for (int i = M - 1;i>=0;i--) {
                if (fa[i][x] != fa[i][y]) {
                    res += dp[i][x] + dp[i][y];
                    x = fa[i][x], y = fa[i][y];
                }
            }
            return res + dp[0][x] + dp[0][y];
        };
        int ans = 0, rt = b[0];
        for (int i = 0;i<s;i++) {
            ans += cal(b[i], b[(i + 1) % s]);
            rt = lca(rt, b[i]);
        }
        ans /= 2;
        if (rt <= n) ans ++;
        cout << ans - s << '\n';
    }
    return;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t; cin >> t;
    while (t--) solve();
    return 0;
}
