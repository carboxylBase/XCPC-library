#define DEBUG 0
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
const db PI = acosl(-1);

void solve() {
    int n, k; cin >> n >> k;
    vector<pll> a(n);
    for (auto &v : a) cin >> v.first >> v.second;

    auto cross = [&](pll a, pll b) -> ll {
        return a.first * b.second - b.first * a.second;
    };

    auto dot = [&](pll a, pll b) -> ll {
        return a.first * b.first + a.second * b.second;
    };

    auto half = [](const pll &p) {
        return p.second < 0 || (p.second == 0 && p.first < 0);
    };

    sort(a.begin(), a.end(), [&](const pll &p, const pll &q) {
        int hp = half(p), hq = half(q);
        if (hp != hq) return hp < hq;
        return cross(p, q) > 0;
    });

    for (int i = 0;i<n;i++) a.push_back(a[i]);

    cout << fixed << setprecision(12);

    db ans = 0;
    for (int i = 0;i<n;i++) {
        ll cr = cross(a[i], a[i + k]);
        ll cd = dot(a[i], a[i + k]);

        db cur = atan2l(cr, cd);
        if (cr < 0 || (cr == 0 && cd > 0 && i + k >= n)) cur += 2 * PI;
        ans = max(ans, cur);
    }

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
