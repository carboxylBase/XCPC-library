#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll INF = 5000000000000000000LL;

void solve() {
    int n;
    ll d;
    cin >> n >> d;
    vector<ll> a(n + 1), b(n + 1);
    for (int i = 1; i <= n; i++) cin >> a[i] >> b[i];
    vector<vector<ll>> dp(n + 1, vector<ll>(4, -INF));
    dp[0][0] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= 3; j++) {
            if (dp[i][j] == -INF) continue;
            for (int k = 0; k <= 3; k++) {
                ll usd = min(b[i + 1], a[i + 1] + j - k);
                if (usd < max(0, j - k)) continue;
                ll res = usd * 2 * d;
                if (usd >= 3 - k) {
                    res += (a[i + 1] + j - k - usd) * d;
                }
                res += k;
                dp[i + 1][k] = max(dp[i + 1][k], dp[i][j] + res);
            }
        }
    }
    cout << *max_element(dp[n].begin(), dp[n].end()) << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
