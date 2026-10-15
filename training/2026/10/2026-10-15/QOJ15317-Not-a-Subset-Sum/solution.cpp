#include <bits/stdc++.h>
using namespace std;
using ll = long long;
void solve() {
    int n; cin >> n;
    int m = 1;
    for (int i = 1; i <= n; i++) m *= 3;
    vector<ll> a(1 << n), b(m);
    for (auto &v : a) cin >> v;
    for (int i = 0; i < (1 << n); i++) {
        int pos = 0, bs = 1;
        for (int j = 0; j < n; j++) {
            if ((i >> j) & 1) pos += bs;
            bs *= 3;
        }
        b[pos] = a[i];
    }
    vector<int> bs(n, 1);
    for (int i = 1; i < n; i++) bs[i] = bs[i-1] * 3;
    for (int i = 0; i < n; i++) {
        int pre = 0;
        for (int j = 0; j <= i; j++) pre += bs[j] * 2;
        auto dfs = [&](auto&& self, int dep, int suf) -> void {
            if (dep == i) {
                for (int j = pre; j >= 0; j--) {
                    int k = ((j + suf) / bs[i]) % 3;
                    if (k < 2) b[j + suf - k * bs[i] + 2 * bs[i]] += b[j + suf];
                    else j -= bs[i] - 1;
                }
                return;
            }
            self(self, dep - 1, suf + bs[dep]);
            self(self, dep - 1, suf);
        };
        dfs(dfs, n - 1, 0);
    }
    ll ans = 0;
    for (auto v : b) ans ^= v;
    cout << ans << '\n';
}
int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    solve();
}
