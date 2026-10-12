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

ll qpow(ll base,ll k,ll mod) {
    if (base == 0) return 0;
    ll res = 1;
    base %= mod; base = (base + mod) % mod;
    k %= (mod - 1); k = (k + mod - 1) % (mod - 1);
    while (k) {
        if (k & 1) {
            res *= base; res %= mod;
        }
        k >>= 1;
        base *= base; base %= mod;
    }
    return res;
}

namespace NiYuan {
    vector<ll> fac, invFac, inv;

    void init(int n) {
        fac.resize(n + 1, 1), invFac.resize(n + 1, 1), inv.resize(n + 1, 1);
        for (int i = 1;i<n+1;i++) {
            fac[i] = fac[i - 1] * i % MOD;
        }
        invFac[n] = qpow(fac[n], MOD - 2, MOD);
        for (int i = n - 1;i>=1;i--) {
            invFac[i] = invFac[i + 1] * (i + 1) % MOD;
        }
        for (int i = 2; i <= n; i++) {
            inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;
        }
    }
};