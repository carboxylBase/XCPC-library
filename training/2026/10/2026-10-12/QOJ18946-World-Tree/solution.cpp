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
const int M = 22;

void solve() {
    int n, q; cin >> n >> q;
    vector<vector<pll>> g(n + 1);
    for (int i = 1;i<n;i++) {
        ll u, v, w; cin >> u >> v >> w;
        g[u].push_back(pll(v, w));
        g[v].push_back(pll(u, w));
    }

    vector<vector<int>> fa(M, vector<int>(n + 1, 0));
    vector<ll> dpF(n + 1, -INF), dp(n + 1, -INF);
    vector<vector<ll>> dp1(M, vector<ll>(n + 1, 0));
    vector<int> dep(n + 1, 0);
    vector<array<pll, 3>> mx(n + 1);
    vector<vector<ll>> dist(M, vector<ll>(n + 1, 0));
    
    auto dfs0 = [&](auto&&self, int u) -> void {
        mx[u][0] = mx[u][1] = mx[u][2] = pll(-INF, 0);
        dp[u] = 0;
        for (auto [v, w] : g[u]) {
            if (v == fa[0][u]) continue;
            fa[0][v] = u;
            dep[v] = dep[u] + 1;
            dist[0][v] = w;
            self(self, v);
            pll z = pll(dp[v] + w, v);
            for (int i = 0;i<3;i++) {
                if (z.first > mx[u][i].first) {
                    swap(z, mx[u][i]);
                }
            }
            dp[u] = max(dp[u], dp[v] + w);
        }
    };
    dep[1] = 1;
    dfs0(dfs0, 1);

    auto dfs1 = [&](auto&&self, int u) -> void {
        vector<ll> pre(g[u].size(), -INF), suf(g[u].size(), -INF);
        for (int i = 0;i<g[u].size();i++) {
            if (i > 0) pre[i] = pre[i - 1];
            if (g[u][i].first == fa[0][u]) continue;
            pre[i] = max(pre[i], dp[g[u][i].first] + g[u][i].second);
        }
        for (int i = int(g[u].size()) - 1;i>=0;i--) {
            if (i + 1 < g[u].size()) suf[i] = suf[i + 1];
            if (g[u][i].first == fa[0][u]) continue;
            suf[i] = max(suf[i], dp[g[u][i].first] + g[u][i].second);
        }
        for (int i = 0;i<g[u].size();i++) {
            if (g[u][i].first == fa[0][u]) continue;
            ll res = 0;
            if (i > 0) res = max(res, pre[i - 1]);
            if (i + 1 < g[u].size()) res = max(res, suf[i + 1]);
            dp1[0][g[u][i].first] = res;
            dpF[g[u][i].first] = max(res, dpF[u]) + g[u][i].second;
            dpF[g[u][i].first] = max(0ll, dpF[g[u][i].first]);
        }
        for (auto [v, w] : g[u]) {
            if (v == fa[0][u]) continue;
            self(self, v);
        }
    };
    dfs1(dfs1, 1);

    for (int i = 1;i<M;i++) {
        for (int j = 1;j<=n;j++) {
            fa[i][j] = fa[i - 1][fa[i - 1][j]];
        }
    }
    for (int i = 1;i<M;i++) {
        for (int j = 1;j<=n;j++) {
            dist[i][j] = dist[i - 1][j] + dist[i - 1][fa[i - 1][j]];
        }
    }
    for (int i = 1;i<M;i++) {
        for (int j = 1;j<=n;j++) {
            dp1[i][j] = max(dp1[i - 1][j], dp1[i - 1][fa[i - 1][j]]);
        }
    }

    auto lca = [&](int x, int y) -> int {
        for (int i = M - 1;i>=0;i--) {
            if (dep[x] > dep[y]) swap(x, y);
            if (dep[fa[i][y]] >= dep[x]) y = fa[i][y];
        }
        if (x == y) return x;
        for (int i = M - 1;i>=0;i--) {
            if (fa[i][x] != fa[i][y]) x = fa[i][x], y = fa[i][y];
        }
        return fa[0][x];
    };

    auto dis = [&](int x, int y) -> ll {
        ll ans = 0;
        for (int i = M - 1;i>=0;i--) {
            if (dep[x] > dep[y]) swap(x, y);
            if (dep[fa[i][y]] >= dep[x]) {
                ans += dist[i][y];
                y = fa[i][y];
            }
        }
        if (x == y) return ans;
        for (int i = M - 1;i>=0;i--) {
            if (fa[i][x] != fa[i][y]) {
                ans += dist[i][x] + dist[i][y];
                x = fa[i][x], y = fa[i][y];
            }
        }
        return ans + dist[0][x] + dist[0][y];
    };

    while (q --) {
        int x, y; cin >> x >> y;

        ll res = dis(x, y);

        if (x == y) {
            cout << 2 * max({0ll, dpF[x], dp[x]}) << '\n';
            continue;
        }

        ll ans = 0;

        int z = lca(x, y);

        if (z == x) {
            ans = max(dp[y], dpF[x]);
        } else if (z == y) {
            ans = max(dp[x], dpF[y]);
        } else {
            ans = max({dp[x], dp[y], dpF[z]});
        }

        for (int i = M - 1;i>=0;i--) {
            if (dep[x] > dep[y]) swap(x, y);
            if (dep[fa[i][y]] > dep[x]) {
                ans = max(ans, dp1[i][y]);
                y = fa[i][y];
            }
        }

        if (dep[x] > dep[y]) swap(x, y);
        if (fa[0][y] == x) {
            for (int i = 0;i<3;i++) {
                if (mx[x][i].second != y) {
                    ans = max(ans, mx[x][i].first);
                    break;
                }
            }
            ans = ans * 2 + res;
            cout << ans << '\n';
            continue;
        }

        if (dep[x] != dep[y]) {
            ans = max(ans, dp1[0][y]);
            y = fa[0][y];
        }

        for (int i = M - 1;i>=0;i--) {
            if (fa[i][x] != fa[i][y]) {
                ans = max({ans, dp1[i][x], dp1[i][y]});
                x = fa[i][x], y = fa[i][y];
            }
        }

        for (int i = 0;i<3;i++) {
            if (mx[z][i].second != x && mx[z][i].second != y) {
                ans = max(ans, mx[z][i].first);
                break;
            }
        }

        cout << ans * 2 + res << '\n';
    }
    return;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
