#include <bits/stdc++.h>
using namespace std;

void solve() {
    int n;
    cin >> n;
    vector<vector<int>> g(n + 1);
    vector<int> fa(n + 1, 0);
    for (int i = 2; i <= n; i++) {
        cin >> fa[i];
        g[fa[i]].push_back(i);
    }
    vector<int> dep(n + 1, 0);
    dep[1] = 1;
    auto dfs = [&](auto&& self, int u) -> void {
        for (auto v : g[u]) {
            dep[v] = dep[u] + 1;
            self(self, v);
        }
    };
    dfs(dfs, 1);

    int lo = 0, hi = n, d;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        cout << "? " << 1 << ' ' << mid << endl;
        int x;
        cin >> x;
        if (x == 1) {
            d = mid;
            hi = mid - 1;
        } else lo = mid + 1;
    }

    vector<int> q;
    for (int i = 1; i <= n; i++) {
        if (dep[i] == d + 1) q.push_back(i);
    }
    while (q.size() > 1) {
        vector<int> dp(n + 1, 0);
        for (auto v : q) dp[v] = 1;
        auto dfs = [&](auto&& self, int u) -> void {
            for (auto v : g[u]) {
                self(self, v);
                dp[u] += dp[v];
            }
        };
        dfs(dfs, 1);

        int mn = 100000000, pos;
        for (int i = 1; i <= n; i++) {
            if (abs(int(q.size()) / 2 - dp[i]) < mn) {
                mn = abs(int(q.size()) / 2 - dp[i]);
                pos = i;
            }
        }
        cout << "? " << pos << ' ' << d + 1 - dep[pos] << endl;
        int x;
        cin >> x;
        if (x == 1) {
            vector<int> nq;
            auto dfs = [&](auto&& self, int u) -> void {
                if (dp[u] == 1 && dep[u] == d + 1) nq.push_back(u);
                for (auto v : g[u]) self(self, v);
            };
            dfs(dfs, pos);
            swap(nq, q);
        } else {
            auto dfs = [&](auto&& self, int u) -> void {
                if (dp[u] == 1 && dep[u] == d + 1) dp[u] = 0;
                for (auto v : g[u]) self(self, v);
            };
            dfs(dfs, pos);
            vector<int> nq;
            auto dfs1 = [&](auto&& self, int u) -> void {
                if (dp[u] == 1 && dep[u] == d + 1) nq.push_back(u);
                for (auto v : g[u]) self(self, v);
            };
            dfs1(dfs1, 1);
            swap(nq, q);
        }
    }
    cout << "! " << q[0] << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t = 1;
    cin >> t;
    while (t--) solve();
    return 0;
}
