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
    ll n, m, P; cin >> n >> m >> P;
    vector<ll> b(n + 1, 0);
    for (int i = 1;i<n+1;i++) cin >> b[i];
    vector<vector<pll>> g(n + 1);
    vector<pll> d;
    for (int i = 1;i<m+1;i++) {
        ll u, v, c; cin >> u >> v >> c;
        g[u].push_back(pll(v, c));
        g[v].push_back(pll(u, c));
        d.push_back(pll(u, (c - b[v] + P) % P));
        d.push_back(pll(v, (c - b[u] + P) % P));
    }
    d.push_back(pll(1, 0));

    sort(d.begin(), d.end());
    d.erase(unique(d.begin(), d.end()), d.end());
    auto idx = [&](pll x) -> int {
        return lower_bound(d.begin(), d.end(), x) - d.begin();
    };

    struct Node {
        int idx, t;
        bool operator<(const Node& A) const {
            return t > A.t;
        }
    };
    vector<ll> dp(d.size(), INF);

    
    priority_queue<Node> nodes;
    Node z;
    z.idx = idx(pll(1, 0));
    z.t = 0;
    dp[z.idx] = z.t;
    nodes.push(z);
    while (!nodes.empty()) {
        Node u = nodes.top();
        nodes.pop();

        if (dp[u.idx] < u.t) continue;

        for (auto [v, w] : g[d[u.idx].first]) {
            pll nxt;
            nxt.first = v;
            nxt.second = (w - b[d[u.idx].first] + P) % P;
            int id = idx(nxt);
            ll dis = (d[u.idx].second - nxt.second + P) % P;
            dis = min(dis, P - dis);
            if (dp[id] > dis + u.t) {
                dp[id] = dis + u.t;
                Node z;
                z.idx = id;
                z.t = dp[id];
                nodes.push(z);
            }
        }
    }

    ll ans = INF;
    for (int i = 0;i<d.size();i++) {
        if (d[i].first == n) {
            ans = min(ans, dp[i]);
        }
    }
    if (ans == INF) ans = -1;
    cout << ans << endl;
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
