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

int half(pll z) {
    return z.second < 0 || (z.second == 0 && z.first < 0);
}

ll cross(pll a, pll b) {
    return a.first * b.second - a.second * b.first;
}

const db PI = acosl(-1);

void solve() {
    ll n, k; cin >> n >> k;
    vector<pll> a(n);
    for (auto &v : a) cin >> v.first >> v.second;
    sort(a.begin(), a.end(), [&](pll a, pll b) {
        if (half(a) != half(b)) return half(a) < half(b);
        return cross(a, b) > 0;
    });

    vector<array<ll, 3>> b;
    for (auto v : a) {
        if (!b.empty() && half(pll(b.back()[0], b.back()[1])) == half(v) && cross(pll(b.back()[0], b.back()[1]), v) == 0) {
            if (abs(b.back()[0]) < abs(v.first) || abs(b.back()[1]) < abs(v.second)) {
                b[b.size() - 1][0] = v.first, b[b.size() - 1][1] = v.second;
            }
            b[b.size() - 1][2] ++;
            continue;
        }
        array<ll, 3> z;
        z[0] = v.first, z[1] = v.second, z[2] = 1;
        b.push_back(z);
    }

    db ans = INF;
    int len = b.size();
    for (int i = 0;i<len;i++) b.push_back(b[i]);
    multiset<ll> mx;
    mx.insert(b[0][0] * b[0][0] + b[0][1] * b[0][1]);
    int l = 0, r = 0;
    int tot = b[0][2];
    for (;l<len;l++) {
        while (r + 1 < b.size() && tot < k) {
            r ++;
            tot += b[r][2];
            mx.insert(b[r][0] * b[r][0] + b[r][1] * b[r][1]);
        }
        if (tot >= k) {
            pll L = {b[l][0], b[l][1]};
            pll R = {b[r][0], b[r][1]};
            ll cr = cross(L, R);
            ll dt = L.first * R.first + L.second * R.second;
            db sita = atan2l((db)cr, (db)dt);
            if (sita < 0) sita += 2 * PI;
            ll r = *mx.rbegin();
            db res = abs(sita) * r / 2;
            ans = min(ans, res);
        }
        tot -= b[l][2];
        mx.erase(mx.find(b[l][0] * b[l][0] + b[l][1] * b[l][1]));
    }
    cout << fixed << setprecision(12) << ans << endl;
    return;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t = 1;
    cin >> t;
    while (t--) solve();
    return 0;
}
