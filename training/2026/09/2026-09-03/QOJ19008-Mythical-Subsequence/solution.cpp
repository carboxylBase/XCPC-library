#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

void solve() {
    int n; cin >> n;
    string s; cin >> s;

    struct Node {
        vector<int> pos;
        int msk;
    };

    vector<vector<Node>> nodes(n);
    map<pii, bool> vis;

    for (int i = 0;i<n;i++) {
        if (s[i] == '0' || s[i] == '9') continue;
        Node z;
        z.msk = (1 << (s[i] - '0'));
        if (!vis.count(pii(z.msk, 1))) {
            z.pos.push_back(i);
            nodes[i].push_back(z);
            vis[pii(z.msk, z.pos.size())] = 1;
        }
        for (int j = i-1;j>=0;j--) {
            for (auto &v : nodes[j]) {
                int sta = 0;
                bool ok = 1;
                sta |= (1 << (s[i] - '0'));
                for (int k = 0;k<9;k++) {
                    if ((v.msk >> k) & 1) {
                        int val = (k * 10 + (s[i] - '0')) % 9;
                        if (val == 0) {
                            ok = 0;
                            break;
                        }
                        sta |= (1 << val);
                    }
                }
                if (!ok) continue;
                if (!vis.count(pii(sta, v.pos.size() + 1))) {
                    Node z;
                    z.pos = v.pos;
                    z.pos.push_back(i);
                    z.msk = sta;
                    nodes[i].push_back(z);
                    vis[pii(sta, v.pos.size() + 1)] = 1;
                }
            }

            if (s[j] == s[i]) {
                break;
            }
        }
    }

    int mx = 0;
    for (int i = 0;i<n;i++) {
        for (auto v : nodes[i]) {
            mx = max(mx, int(v.pos.size()));
        }
    }
    for (int i = 0;i<n;i++) {
        for (auto v : nodes[i]) {
            if (mx == (int)v.pos.size()) {
                cout << mx << ' ';
                for (auto z : v.pos) {
                    cout << z + 1 << ' ';
                }
                cout << '\n';
                return;
            }
        }
    }
    if (mx == 0) {
        cout << mx << '\n';
        return;
    }
    return;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t = 1;
    cin >> t;

    while (t--) {
        solve();
    }

    return 0;
}
