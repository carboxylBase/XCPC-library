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

vector<int> q = {1, 2, 5, 10, 20, 50, 100, 200, 500};
vector<int> b = {4 ,40 ,44 ,400 ,404 ,440 ,444 ,1200 ,1204 ,1240 ,1244 ,2000 ,2004 ,2040 ,2044 ,2400 ,2404 ,2440 ,2444 ,3200 ,3204 ,3240 ,3244 ,4000 ,4004 ,4040 ,4044 ,4400 ,4404 ,4440 ,4444 ,5200 ,5204 ,5240 ,5244 ,6000 ,6004 ,6040 ,6044 ,6400 ,6404 ,6440 ,6444 ,7200 ,7204 ,7240 ,7244 ,8000 ,8004 ,8040 ,8044 ,8400 ,8404 ,8440 ,8444 ,9200 ,9204 ,9240 ,9244 ,10000 };

vector<int> chk(int x) {
    vector<int> ans;
    auto dfs = [&](auto&&self, vector<bool> vis, vector<int> cnt, int dep) -> void {
        if (dep == q.size()) return;
        if (!ans.empty()) return;

        for (int i = 0;i*q[dep]<=x;i++) {
            bool bad = 0;
            for (auto v : q) {
                if (v == q[dep]) continue;
                if (v == i * q[dep]) {
                    bad = 1;
                    break;
                }
            }
            if (bad) break;
            vector<bool> nvis = vis;
            for (int j = x;j>=0;j--) {
                for (int k = 1;k<=i;k++) {
                    if (nvis[j] && j + k * q[dep] <= x) {
                        nvis[j + k * q[dep]] = 1;
                    }
                }
            }
            if (nvis[x / 2]) continue;
            if (nvis[x]) {
                ans = cnt;
                return;
            } else {
                vector<int> ncnt = cnt;
                ncnt.push_back(i);
                self(self, nvis, ncnt, dep + 1);
            }
        }
    };
    vector<bool> vis(x + 1, 0);
    vis[0] = 1;
    dfs(dfs, vis, vector<int>(0), 0);
    return ans;
}

void solve() {
    for (auto v : b) {
        vector<int> ans = chk(v);
        while (ans.size() < q.size()) ans.push_back(0);
        cout << "ans[v] = {";
        for (int i = 0;i<ans.size();i++) {
            cout << ans[i];
            if (i != ans.size() - 1) {
                cout << ", ";
            }
        }
        cout << "};";
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
