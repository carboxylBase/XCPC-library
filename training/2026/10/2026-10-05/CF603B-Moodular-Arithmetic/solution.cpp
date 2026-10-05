#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const int N = 2000000;
const ll MOD = 1000000007;
int fa[N];

int f(int x) {
    if (fa[x] == x) return x;
    return fa[x] = f(fa[x]);
}

void merge(int x, int y) {
    x = f(x);
    y = f(y);
    fa[x] = y;
}

void solve() {
    int p, k;
    cin >> p >> k;
    if (k == 0) {
        ll ans = 1;
        for (int i = 1; i < p; i++) ans = ans * p % MOD;
        cout << ans << '\n';
        return;
    }

    for (int i = 0; i < p; i++) fa[i] = i;
    for (int i = 1; i < p; i++) merge(i, 1LL * i * k % p);

    ll ans = 1;
    for (int i = 1; i < p; i++) {
        if (f(i) == i) ans = ans * p % MOD;
    }
    if (k == 1) ans = ans * p % MOD;
    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
