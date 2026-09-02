#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pii = pair<int, int>;

const ll MOD = 1e9 + 7;

namespace ComNum {
    const int MAXN = 200005;
    ll fac[MAXN], inv[MAXN];

    ll qpow(ll base, ll k, ll mod) {
        if (base == 0) return 0;
        ll res = 1;
        base %= mod;
        base = (base + mod) % mod;
        k %= mod - 1;
        k = (k + mod - 1) % (mod - 1);
        while (k) {
            if (k & 1) {
                res *= base;
                res %= mod;
            }
            k >>= 1;
            base *= base;
            base %= mod;
        }
        return res;
    }

    void init() {
        fac[0] = 1;
        for (int i = 1; i < MAXN; i++) {
            fac[i] = fac[i - 1] * i % MOD;
        }

        vector<ll> pre(MAXN, 1), suf(MAXN, 1);
        for (int i = 0; i < MAXN; i++) {
            if (i - 1 >= 0) pre[i] = pre[i - 1];
            pre[i] = pre[i] * fac[i] % MOD;
        }
        for (int i = MAXN - 1; i >= 0; i--) {
            if (i + 1 < MAXN) suf[i] = suf[i + 1];
            suf[i] = suf[i] * fac[i] % MOD;
        }

        ll base = qpow(pre.back(), MOD - 2, MOD);
        for (int i = 0; i < MAXN; i++) {
            ll numerator = 1;
            if (i - 1 >= 0) numerator = pre[i - 1] * numerator % MOD;
            if (i + 1 < MAXN) numerator = suf[i + 1] * numerator % MOD;
            inv[i] = numerator * base % MOD;
        }
    }

    ll C(ll n, ll m) {
        if (m == 0) return 1;
        if (n >= 0) {
            if (n < m) return 0;
            return fac[n] * inv[n - m] % MOD * inv[m] % MOD;
        }

        ll a = -n;
        ll nn = a + m - 1;
        ll res = fac[nn] * inv[nn - m] % MOD * inv[m] % MOD;
        if (m & 1) res = (MOD - res) % MOD;
        return res;
    }
}

void solve() {
    int h, w, n;
    cin >> h >> w >> n;
    vector<pii> points(n);
    for (auto &point : points) {
        cin >> point.first >> point.second;
    }

    sort(points.begin(), points.end());
    if (points.back() != pii(h, w)) {
        points.push_back(pii(h, w));
        n++;
    }

    vector<ll> dp(n, 0);
    for (int i = 0; i < n; i++) {
        dp[i] = ComNum::C(points[i].first + points[i].second - 2,
                          points[i].first - 1);
        for (int j = 0; j < i; j++) {
            if (points[j].first <= points[i].first &&
                points[j].second <= points[i].second) {
                dp[i] -= dp[j] * ComNum::C(
                    points[i].first - points[j].first +
                        points[i].second - points[j].second,
                    points[i].first - points[j].first) % MOD;
                dp[i] = (dp[i] + MOD) % MOD;
            }
        }
    }

    cout << dp.back() << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ComNum::init();
    solve();
    return 0;
}
