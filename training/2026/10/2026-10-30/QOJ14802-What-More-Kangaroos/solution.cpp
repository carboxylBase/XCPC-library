#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using pll = pair<ll,ll>;

void solve() {
    int n; cin >> n;
    vector<pll> a;
    for (int i = 0; i < n; i++) {
        pll z; cin >> z.first >> z.second;
        int x; cin >> x;
        if (z.first != 0 || z.second != 0) a.push_back(z);
    }
    n = a.size();

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
    for (int i = 0; i < n; i++) a.push_back(a[i]);

    int ans = 0;
    for (int i = 0, j = 0; i < n; i++) {
        if (j < i) j = i;
        while (j + 1 - i < n) {
            ll val_r = cross(a[i], a[j + 1]);
            ll val_d = dot(a[j + 1], a[i]);
            if (val_r == 0 && val_d > 0 && j + 1 >= n) break;
            if (val_r == 0 && val_d < 0) break;
            if (val_r < 0) break;
            j++;
        }
        ans = max(ans, j - i + 1);
    }
    cout << ans << endl;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t = 1;
    cin >> t;
    while (t--) solve();
    return 0;
}
