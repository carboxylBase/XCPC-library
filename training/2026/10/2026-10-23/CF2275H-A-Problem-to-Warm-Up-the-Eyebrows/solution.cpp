#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 1e9 + 7;

void solve() {
    int n, m; cin >> n >> m;
    vector<vector<ll>> a(n + 1, vector<ll>(m + 1, 0));
    for (int i = 1;i<n+1;i++) {
        for (int j = 1;j<m+1;j++) {
            cin >> a[i][j];
            a[i][j] = (a[i][j] % MOD + MOD) % MOD;
        }
    }
    ll T = 0, D = 0;
    vector<vector<ll>> b(n + 2, vector<ll>(m + 2, 0));
    for (int i = 1;i<n+1;i++) {
        for (int j = 1;j<m+1;j++) {
            if ((i & 1) && (j & 1)) b[i][j] = a[i][j];
            b[i][j] += b[i - 1][j] + b[i][j - 1] - b[i - 1][j - 1];
            b[i][j] = (b[i][j] % MOD + MOD) % MOD;
        }
    }
    for (int i = 1;i<n+1;i++) {
        for (int j = 1;j<m+1;j++) {
            if (i % 2 == n % 2 && j % 2 == m % 2)
                T = (T + a[i][j] * b[i][j] % MOD) % MOD;
        }
    }
    for (int i = 1;i<n+1;i++) {
        for (int j = m;j>=1;j--) {
            b[i][j] = 0;
            if ((i & 1) && j % 2 == m % 2) b[i][j] = a[i][j];
            b[i][j] += b[i - 1][j] + b[i][j + 1] - b[i - 1][j + 1];
            b[i][j] = (b[i][j] % MOD + MOD) % MOD;
        }
    }
    for (int i = 1;i<n+1;i++) {
        for (int j = m;j>=1;j--) {
            if (i % 2 == n % 2 && (j & 1))
                T = (T + a[i][j] * b[i - 1][j + 1] % MOD) % MOD;
        }
    }
    for (int i = 1;i<n+1;i++) {
        for (int j = 1;j<m+1;j++) {
            if ((i & 1) && i % 2 == n % 2) {
                if ((j & 1) && j % 2 == m % 2)
                    D = (D + a[i][j] * a[i][j] % MOD) % MOD;
            }
        }
    }
    ll ans = (2 * T - D) % MOD;
    if ((n + m) & 1) ans = (-ans) % MOD;
    ans = (ans + MOD) % MOD;
    cout << ans << endl;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t = 1;
    cin >> t;
    while (t--) solve();
    return 0;
}
