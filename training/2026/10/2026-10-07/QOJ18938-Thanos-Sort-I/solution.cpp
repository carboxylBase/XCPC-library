#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 998244353;

ll qpow(ll base, ll k) {
    ll res = 1;
    while (k) {
        if (k & 1) res = res * base % MOD;
        base = base * base % MOD;
        k >>= 1;
    }
    return res;
}

void solve() {
    int n;
    cin >> n;
    vector<ll> fac(n + 1, 1);
    for (int i = 1; i <= n; i++) fac[i] = fac[i - 1] * i % MOD;
    auto cal = [&](auto&& self, int x) -> ll {
        if (x == 1) return 1;
        ll res = (self(self, (x + 1) / 2)
                  + qpow(fac[x - 1], MOD - 2)
                  - (x + 1) / 2 * qpow(fac[x], MOD - 2) % MOD) % MOD;
        return (res + MOD) % MOD;
    };
    cout << cal(cal, n) << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
