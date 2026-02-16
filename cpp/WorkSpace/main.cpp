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

void solve() {
    struct F {
        ll a, b, c;
    };
    int n; cin >> n;
    vector<F> a(n);
    for (auto &v : a) {
        cin >> v.a >> v.b >> v.c;
    }

    auto chk = [&](F& L, F& R) -> bool {
        if (L.a == R.a) {
            if (L.b == R.b && L.c != R.c) {
                return 0;
            }
        } else {
            if ((L.b - R.b) * (L.b - R.b) - 4 * (L.a - R.a) * (L.c - R.c) < 0) {
                return 0;
            }
        }
        return 1;
    };

    auto cal = [&](F& L) -> db {
        return -db(L.b * L.b) / 4 / L.a + L.c;
    };

    auto baohan = [&](F& L, F& R) -> bool {
        if (L.a * R.a < 0) return 0;
        if (chk(L, R)) return 0;
        db l = cal(L), r = cal(R);
        if (L.a > 0) {
            if (l < r) return 1;
        } else {
            if (l > r) return 1;
        }
        return 0;
    };

    vector<vector<int>> g(n), e(n);
    for (int i = 0;i<n;i++) {
        for (int j = 0;j<n;j++) {
            if (i == j) continue;
            if (baohan(a[i], a[j])) {
                g[i].push_back(j);
                e[j].push_back(i);
            }
        }
    }

    queue<int> q;
    vector<int> in(n, 0);
    vector<int> dp1(n, 0);
    for (int i = 0;i<n;i++) {
        in[i] = g[i].size();
        if (in[i] == 0) {
            q.push(i);
        }
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto v : g[u]) {
            dp1[u] = max(dp1[u], dp1[v]);
        }
        dp1[u] ++;
        for (auto v : e[u]) {
            in[v] --;
            if (in[v] == 0) q.push(v);
        }
    }

    vector<int> dp2(n, 0);
    for (int i = 0;i<n;i++) {
        for (int j = 0;j<n;j++) {
            if (i == j) continue;
            if (a[i].a * a[j].a > 0) continue;
            if (chk(a[i], a[j])) continue;
            dp2[i] = max(dp2[i], dp1[j]);
        }
    }

    vector<int> dp3(n, 0);
    for (int i = 0;i<n;i++) {
        in[i] = e[i].size();
        if (in[i] == 0) q.push(i);
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto v : e[u]) {
            dp3[u] = max(dp3[u], dp3[v] + 1);
        }
        dp3[u] = max(dp3[u], dp2[u]);
        for (auto v : g[u]) {
            in[v] --;
            if (in[v] == 0) q.push(v);
        }
    }

    // for (int i = 0;i<n;i++) {
    //     cout << i << ' ' << dp1[i] << ' ' << dp2[i] << ' ' << dp3[i] << endl;         
    // }
    // cout << endl;

    for (int i = 0;i<n;i++) {
        cout << dp3[i] + dp1[i] << ' ';
    }
    cout << endl;
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
    cin >> t;

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
