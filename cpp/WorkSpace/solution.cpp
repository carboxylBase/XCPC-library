#define DEBUG 0
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

vector<array<int, 6>> val = {
    {0, 0, 1, 0, 1, 1}, 
    {0, 1, 0, 2, 1, 2}};

void solve() {
    int n; cin >> n;
    vector<array<int, 6>> ans;
    
    if (3 * n % 6 == 0) {
        int cur = 1;
        while (cur <= 3 * n) {
            int pos = 1;
            while (pos <= 3 * n) {
                for (auto u : val) {
                    for (int i = 0;i<6;i++) {
                        if (i & 1) {
                            u[i] += pos;
                        } else {
                            u[i] += cur;
                        }
                    }
                    ans.push_back(u);
                }
                pos += 3;
            }
            cur += 2;
        }
    } else {
        int cur = 1;
        while (cur <= 3 * n) {
            int pos = 1;
            while (pos <= 3 * n) {
                for (auto u : val) {
                    for (int i = 0;i<6;i++) {
                        if (i & 1) {
                            u[i] += pos;
                        } else {
                            u[i] += cur;
                        }
                    }
                    ans.push_back(u);
                }
                pos += 3;
            }
            cur += 2;
        }
        for (int i = 0;i<4*n;i++) {
            ans.pop_back();
        }
        array<int, 6> delta = {0, 0, 1, 0, 2, 1};
        int x = 3 * n - 2, y = 1;
        for (int i = 1;i<3*n;i++) {
            array<int, 6> z = delta;
            for (int j = 0;j<6;j++) {
                if (j & 1) {
                    z[j] += y;
                } else {
                    z[j] += x;
                }
            }
            ans.push_back(z);
            y ++;
        }
    }

    {
        vector<array<int, 6>> res;
        for (auto v : ans) {
            bool ok = 1;
            for (auto u : v) {
                if (u > 3 * n || u < 1) {
                    ok = 0;
                }
            }
            if (ok) res.push_back(v);
        }
        ans.swap(res);
    }

    cout << ans.size() << endl;
    for (auto u : ans) {
        for (auto v : u) {
            cout << v << ' ';
        }
        cout << endl;
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
    cin >> t;

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
