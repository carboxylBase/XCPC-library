#include <bits/stdc++.h>
using namespace std;
using ll = long long;

struct LinearBasis {
    ll basis[63];

    void insert(ll x) {
        if (!x) return;
        for (int i = 60; i >= 0; i--) {
            if ((1ll << i) & x) {
                if (basis[i]) x ^= basis[i];
                else {
                    basis[i] = x;
                    break;
                }
            }
        }
    }

    void preprocess() {
        for (int i = 60; i >= 0; i--) {
            if (!basis[i]) continue;
            for (int j = i - 1; j >= 0; j--) {
                if (basis[j] && ((1ll << j) & basis[i])) {
                    basis[i] ^= basis[j];
                }
            }
        }
    }

    void init() {
        memset(basis, 0, sizeof basis);
    }
} solver, s1;

void solve() {
    int n;
    cin >> n;
    vector<ll> a(n);
    for (auto &v : a) cin >> v;

    solver.init();
    for (auto v : a) solver.insert(v);
    solver.preprocess();

    auto chk = [&](ll val) -> bool {
        ll res = 0;
        for (auto v : a) res ^= v;
        res &= val;
        if (res == val) return true;
        if (res == 0) {
            s1.init();
            for (int i = 0; i < 60; i++) {
                if (solver.basis[i]) {
                    ll res = solver.basis[i];
                    res &= val;
                    s1.insert(res);
                }
            }
            s1.preprocess();
            for (int i = 59; i >= 0; i--) {
                if ((val >> i) & 1) val ^= s1.basis[i];
            }
            return val == 0;
        }
        return false;
    };

    ll ans = 0;
    for (int i = 59; i >= 0; i--) {
        if (chk(ans + (1ll << i))) ans += (1ll << i);
    }
    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    cin >> t;
    while (t--) solve();
    return 0;
}
