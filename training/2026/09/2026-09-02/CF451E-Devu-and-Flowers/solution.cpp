#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pll = pair<ll, ll>;

const ll MOD = 1e9 + 7;

namespace ComNum {
    const int MAXN = 25;
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
        } else {
            n = -n;
            ll res = 1;
            for (int i = 0; i < n - 1; i++) {
                res = res * ((n + m - 1 - i) % MOD) % MOD;
            }
            res = res * inv[n - 1] % MOD;
            if (m & 1) res = (MOD - res) % MOD;
            return res;
        }
    }
}

void solve() {
    ll n, s;
    cin >> n >> s;
    vector<ll> f(n + 1, 0);
    for (int i = 1; i <= n; i++) cin >> f[i];

    vector<pll> terms(1, pll(1, 0));
    for (int i = 1; i <= n; i++) {
        vector<pll> next;
        next.reserve(terms.size() * 2);
        for (auto term : terms) {
            next.push_back(pll(term.first, term.second));
            next.push_back(pll(-term.first, term.second + f[i] + 1));
        }
        terms = move(next);
    }

    ll ans = 0;
    for (auto term : terms) {
        ll t = s - term.second;
        if (t < 0) continue;
        ll res = ComNum::C(-n, t);
        if (t & 1) res = (MOD - res) % MOD;
        res = res * term.first % MOD;
        res = (res + MOD) % MOD;
        ans = (ans + res) % MOD;
    }

    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ComNum::init();
    solve();
    return 0;
}

