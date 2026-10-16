#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll INF = 5000000000000000000LL;

void solve() {
    ll n, C;
    cin >> n >> C;
    vector<ll> a(n + 2, 0);
    for (int i = 1; i <= n; i++) cin >> a[i];

    vector<vector<ll>> dp(1 << 9, vector<ll>(1 << 9, INF));
    auto get_high = [&](int x) -> int { return x >> 9; };
    int z = (1 << 9) - 1;
    auto get_low = [&](int x) -> int { return z & x; };

    auto upd = [&](int pos, ll val) -> void {
        int h = get_high(pos), l = get_low(pos);
        for (int i = 0; i < (1 << 9); i++) {
            dp[h][i] = min(dp[h][i], (i ^ l) + val);
        }
    };
    upd(0, 0);

    for (int i = 1; i <= n + 1; i++) {
        int h = get_high(a[i]);
        int l = get_low(a[i]);
        ll val = INF;
        for (int j = 0; j < (1 << 9); j++) {
            val = min(val, dp[j][l] + ((h ^ j) << 9));
        }
        val += (i - 1) * C;
        upd(a[i], val - i * C);
        if (i == n + 1) cout << val << '\n';
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    cin >> t;
    while (t--) solve();
    return 0;
}
