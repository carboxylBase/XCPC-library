#define DEBUG 1
#define FUCK cout << "fuck" << endl;
#if DEBUG
    #include "all.hpp"
#else
    #include <bits/stdc++.h>
#endif

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

    auto dfs = [&](auto&& self, int u) -> void {
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
        if (fa[0][i] <= n) {
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
        int u, v; cin >> u >> v;
        int ans = 0;
        for (int i = M - 1; i >= 0; i--) {
            if (dep[u] != dep[v]) {
                if (dep[v] > dep[u]) swap(u, v);
                if (dep[fa[i][u]] >= dep[v]) {
                    ans += dp[i][u];
                    u = fa[i][u];
                }
            }
        }
        if (u == v) {
            ans++;
            cout << ans << '\n';
            continue;
        }
        for (int i = M - 1; i >= 0; i--) {
            if (fa[i][u] != fa[i][v]) {
                ans += dp[i][u] + dp[i][v];
                u = fa[i][u], v = fa[i][v];
            }
        }
        ans += 2;
        if (fa[0][u] <= n) ans++;
        cout << ans << '\n';
    }
    return;
}

signed main() {
#if DEBUG
    freopen("input.txt", "r", stdin);
    auto start_time = chrono::steady_clock::now();
#else
    ios::sync_with_stdio(false);
#endif
    cin.tie(nullptr);

    int t = 1;
    // cin >> t;

    while (t--) {
        solve();
    }

#if DEBUG
    auto end_time = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cerr << "Time: " << diff.count() << " ms" << endl;

#endif

    return 0;
}
