#include <bits/stdc++.h>
using namespace std;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<vector<int>> a(n + 1);
    for (int i = 1; i <= m; i++) {
        int k;
        cin >> k;
        while (k--) {
            int x;
            cin >> x;
            a[x].push_back(i);
        }
    }
    const int dummy = m + 1;
    vector<vector<int>> g(m + 2);
    for (int i = 1; i <= n; i++) {
        if (a[i].size() == 1) {
            g[a[i][0]].push_back(dummy);
            g[dummy].push_back(a[i][0]);
        } else if (a[i].size() == 2) {
            g[a[i][0]].push_back(a[i][1]);
            g[a[i][1]].push_back(a[i][0]);
        }
    }
    int ans = 0;
    for (int i = 1; i <= m; i++) ans += min(2, (int)g[i].size());
    vector<bool> vis(m + 2, false);
    vector<int> s;
    for (int i = 1; i <= m; i++) {
        if (vis[i]) continue;
        bool had = false, allTwo = true;
        int cnt = 0;
        s.clear();
        s.push_back(i);
        vis[i] = true;
        while (!s.empty()) {
            int u = s.back();
            s.pop_back();
            cnt++;
            if (u == dummy) had = true;
            if (g[u].size() != 2) allTwo = false;
            for (int v : g[u]) {
                if (!vis[v]) {
                    vis[v] = true;
                    s.push_back(v);
                }
            }
        }
        if (!had && allTwo && (cnt & 1)) ans--;
    }
    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
