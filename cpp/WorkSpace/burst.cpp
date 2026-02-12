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

ll dp[9][82][1<<9];
ll tot[1 << 9][1 << 9];
pii lst[9][82][1<<9];
void solve() {
    int T; cin >> T;
    struct Node {
        int n, idx, k;
    };
    vector<Node> nodes[10];
    vector<ll> val(T);
    vector<vector<string>> ans(T);
    for (int i = 0;i<T;i++) {
        int n, m; cin >> n >> m;
        int k; cin >> k;
        Node z;
        z.k = k;
        z.n = n, z.idx = i;
        nodes[m].push_back(z);
    }
    for (int i = 1;i<=9;i++) {
        for (int o = 0;o<9;o++) {
            for (int j = 0;j<=9*i;j++) {
                for (int k = 0;k<(1<<i);k++) {
                    dp[o][j][k] = -INF;
                }
            }
        }
        for (int p = 0;p<(1<<i);p++) {
            for (int q = 0;q<(1<<i);q++) {
                int cnt = 0;
                for (int x = 0;x<i;x++) {
                    if ((q >> x) & 1) {
                        if (x - 1 >= 0) {
                            if (!((p >> (x - 1)) & 1)) {
                                cnt ++;
                            }
                        }
                        if (x + 1 < i) {
                            if (!((p >> (x + 1)) & 1)) {
                                cnt ++;
                            }
                        }
                        if (!((p >> x) & 1)) {
                            cnt ++;
                        }
                    } else {
                        if (x - 1 >= 0) {
                            if ((q >> (x - 1)) & 1) {
                                cnt ++;
                            }
                            if ((p >> (x - 1)) & 1) {
                                cnt ++;
                            }
                        }
                        if (x + 1 < i) {
                            if ((q >> (x + 1)) & 1) {
                                cnt ++;
                            }
                            if ((p >> (x + 1)) & 1) {
                                cnt ++;
                            }
                        }
                        if ((p >> x) & 1) {
                            cnt ++;
                        }
                    }
                }
                tot[p][q] = cnt;
            }
        }
        for (int j = 0;j<(1<<i);j++) {
            int cnt = 0;
            for (int k = 0;k<i;k++) {
                if ((j >> k) & 1) continue;
                if (k - 1 >= 0) {
                    if ((j >> (k - 1)) & 1) {
                        cnt ++;
                    }
                }
                if (k + 1 < i) {
                    if ((j >> (k + 1)) & 1) {
                        cnt ++;
                    }
                }
            }
            dp[0][__popcount(j)][j] = cnt;
        }
        for (int j = 1;j<9;j++) {
            for (int k = 0;k<=j*i;k++) {
                for (int q = 0;q<(1<<i);q++) {
                    int d = __popcount(q);
                    for (int p = 0;p<(1<<i);p++) {
                        int cnt = tot[p][q];
                        if (dp[j - 1][k][p] + cnt > dp[j][k + d][q]) {
                            dp[j][k + d][q] = dp[j - 1][k][p] + cnt;
                            lst[j][k + d][q] = pii(k, p);
                        }
                    }
                }
            }
        }

        for (auto v : nodes[i]) {
            ll mx = -INF;
            pii sta;
            for (int j = 0;j<(1<<i);j++) {
                if (mx < dp[v.n - 1][v.k][j]) {
                    mx = dp[v.n - 1][v.k][j];
                    sta = pii(v.k, j);
                }
            }
            val[v.idx] = mx;
            for (int j = v.n - 1;j>=0;j--) {
                string z;
                for (int k = 0;k<i;k++) {
                    if ((sta.second >> k) & 1) {
                        z.push_back('*');
                    } else {
                        z.push_back('.');
                    }
                }
                ans[v.idx].push_back(z);
                sta = lst[j][sta.first][sta.second];
            }
        }
    }

    for (int i = 0;i<T;i++) {
        cout << val[i] << '\n';
        for (auto v : ans[i]) {
            cout << v << '\n';
        }
    }
    return;
}

signed main() {
#if DEBUG
    freopen("input.txt", "r", stdin);
    auto start_time = chrono::steady_clock::now();
#else
    ios::sync_with_stdio(false);
#endif
    cin.tie(nullptr);

    int t = 1;
    // cin >> t;

    while (t--) {
        solve();
    }

#if DEBUG
    auto end_time = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cerr << "Time: " << diff.count() << " ms" << endl;

#endif

    return 0;
}
