#include <bits/stdc++.h>

using namespace std;
using ll = long long;

const ll MOD = 998244353;

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
        if (n < m || m < 0) return 0;
        return fac[n] * inv[n - m] % MOD * inv[m] % MOD;
    }
}

void solve() {
    ll n, k;
    cin >> n >> k;

    if (k == 0) {
        cout << ComNum::fac[n] << '\n';
        return;
    }
    if (k >= n) {
        cout << 0 << '\n';
        return;
    }

    ll m = n - k;
    ll ans = 0;

    for (int i = 0; i <= m; i++) {
        ll res = ComNum::C(m, i) * ComNum::qpow(i, n, MOD) % MOD;
        if ((m - i) & 1) {
            res = (MOD - res) % MOD;
        }
        ans = (ans + res) % MOD;
    }

    ans = ans * ComNum::C(n, m) % MOD;
    cout << ans * 2 % MOD << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ComNum::init();
    solve();
    return 0;
}

