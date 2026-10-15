#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;
const int INF = 100000000;
int dp[2][30][30];

void solve() {
    int n, k; cin >> n >> k;
    vector<int> p(n + 1, 0);
    for (int i = 1; i <= n; i++) cin >> p[i];
    vector<bool> vis(n + 1, 0);
    vector<int> a;
    for (int i = 1; i <= n; i++) {
        if (vis[i]) continue;
        int c = 0, pos = i;
        while (!vis[pos]) {
            vis[pos] = 1;
            c++;
            pos = p[pos];
        }
        a.push_back(c);
    }
    sort(a.begin(), a.end());
    if (a.back() > k + 1) {
        cout << "Bob\n";
        return;
    }
    if (a.back() <= k) {
        reverse(a.begin(), a.end());
        int f = 0;
        for (int i = 0; i <= k; i++) {
            for (int j = 0; j <= k; j++) dp[0][i][j] = INF;
        }
        dp[0][0][0] = 0;
        for (int i = 0; i < (int)a.size(); i++) {
            for (int j = 0; j <= k; j++) {
                for (int x = 0; x <= k; x++) dp[f ^ 1][j][x] = INF;
            }
            for (int j = k; j >= 0; j--) {
                for (int x = k; x >= 0; x--) {
                    if (i - j < k) {
                        if (j + 1 <= k && x + a[i] <= k)
                            dp[f ^ 1][j + 1][x + a[i]] = min(dp[f ^ 1][j + 1][x + a[i]], dp[f][j][x]);
                        dp[f ^ 1][j][x] = min(dp[f ^ 1][j][x], dp[f][j][x] + a[i] - 1);
                    } else {
                        dp[f ^ 1][j][x] = min(dp[f ^ 1][j][x], dp[f][j][x]);
                        if (j + 1 <= k && x + a[i] <= k)
                            dp[f ^ 1][j + 1][x + a[i]] = min(dp[f ^ 1][j + 1][x + a[i]], dp[f][j][x]);
                    }
                }
            }
            f ^= 1;
        }
        for (int i = 1; i <= k; i++) {
            for (int j = 0; j <= k; j++) {
                if (dp[f][i][j] <= k) {
                    cout << "Alice\n";
                    return;
                }
            }
        }
        cout << "Bob\n";
        return;
    }
    int c = 0;
    for (auto v : a) if (v == k + 1) c++;
    if (k >= 2 && c >= 2) {
        cout << "Bob\n";
        return;
    }
    if (k >= 2) {
        int one = 0, tot = 0;
        for (auto v : a) {
            if (v == 1) one++;
            else if (v < k + 1) tot += v;
        }
        if (tot > k) {
            cout << "Bob\n";
            return;
        }
        vector<int> dp1(n + 1, 0);
        for (int i = 1; i <= k; i++) dp1[i] = 1;
        for (int i = k + 1; i <= n; i++) {
            for (int j = 1; j <= k; j++) {
                if (dp1[i - j] == 0) dp1[i] = 1;
            }
        }
        for (int i = 0; i <= k - tot && i <= one; i++) {
            if (tot == 0 && i == 0) continue;
            if (dp1[one - i] == 0) {
                cout << "Alice\n";
                return;
            }
        }
        cout << "Bob\n";
        return;
    }
    c = 0;
    for (auto v : a) if (v == 1) c++;
    cout << (c & 1 ? "Alice" : "Bob") << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t; cin >> t;
    while (t--) solve();
    return 0;
}
