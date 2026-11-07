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
vector<int> b = {4 ,40 ,44 ,400 ,404 ,440 ,444 ,1200 ,1204 ,1240 ,1244 ,2000 ,2004 ,2040 ,2044 ,2400 ,2404 ,2440 ,2444 ,3200 ,3204 ,3240 ,3244 ,4000 ,4004 ,4040 ,4044 ,4400 ,4404 ,4440 ,4444 ,5200 ,5204 ,5240 ,5244 ,6000 ,6004 ,6040 ,6044 ,6400 ,6404 ,6440 ,6444 ,7200 ,7204 ,7240 ,7244 ,8000 ,8004 ,8040 ,8044 ,8400 ,8404 ,8440 ,8444 ,9200 ,9204 ,9240 ,9244 ,10000 };
map<int, vector<int>> ans;
vector<int> q = {1, 2, 5, 10, 20, 50, 100, 200, 500};

void solve() {
    int n; cin >> n;

    if (n % 2 == 1) {
        cout << n << endl;
        for (int i = 1;i<=n;i++) {
            cout << 1 << ' ';
        }
        cout << endl;
        return;
    }

    for (auto v : q) {
        if (n % v < (n / 2) % v) {
            vector<int> ans;
            for (int i = 0;i<n%v;i++) {
                ans.push_back(1);
            }
            for (int i = 0;i<n/v;i++) {
                ans.push_back(v);
            }
            cout << ans.size() << endl;
            for (auto v : ans) {
                cout << v << ' ';
            }
            cout << endl;
            return;
        }
    }

    ans[9200] = {0, 0, 0, 0, 0, 0, 0, 1, 18};
    ans[9204] = {0, 2, 0, 0, 0, 0, 0, 1, 18};
    ans[9240] = {0, 0, 0, 0, 2, 0, 0, 1, 18};
    ans[9244] = {0, 2, 0, 0, 2, 0, 0, 1, 18};

    ans[7200] = {0, 0, 0, 0, 0, 0, 0, 1, 14};
    ans[7204] = {0, 2, 0, 0, 0, 0, 0, 1, 14};
    ans[7240] = {0, 0, 0, 0, 2, 0, 0, 1, 14};
    ans[7244] = {0, 2, 0, 0, 2, 0, 0, 1, 14};

    ans[5200] = {0, 0, 0, 0, 0, 0, 0, 1, 10};
    ans[5204] = {0, 2, 0, 0, 0, 0, 0, 1, 10};
    ans[5240] = {0, 0, 0, 0, 2, 0, 0, 1, 10};
    ans[5244] = {0, 2, 0, 0, 2, 0, 0, 1, 10};

    ans[3200] = {0, 0, 0, 0, 0, 0, 0, 1, 6};
    ans[3204] = {0, 2, 0, 0, 0, 0, 0, 1, 6};
    ans[3240] = {0, 0, 0, 0, 2, 0, 0, 1, 6};
    ans[3244] = {0, 2, 0, 0, 2, 0, 0, 1, 6};

    ans[1200] = {0, 0, 0, 0, 0, 0, 0, 1, 2};
    ans[1204] = {0, 2, 0, 0, 0, 0, 0, 1, 2};
    ans[1240] = {0, 0, 0, 0, 2, 0, 0, 1, 2};
    ans[1244] = {0, 2, 0, 0, 2, 0, 0, 1, 2};

    if (ans.count(n)) {
        vector<int> res;
        for (int i = 0;i<ans[n].size();i++) {
            for (int j = 0;j<ans[n][i];j++) {
                res.push_back(q[i]);
            }
        }
        cout << res.size() << endl;
        for (auto v : res) {
            cout << v << ' ';
        }
        cout << endl;
    } else {
        cout << "splittable" << endl;
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
