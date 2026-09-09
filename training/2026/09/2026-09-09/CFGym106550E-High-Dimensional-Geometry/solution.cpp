#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

void solve() {
    int n; cin >> n;
    ll tot = 0;
    for (int i = 1; i < (1 << n); i++) {
        vector<int> q[2];
        for (int j = 0; j < (1 << n); j++) {
            int c = 0;
            for (int k = 0; k < n; k++) {
                if ((i >> k) & 1) {
                    if ((j >> k) & 1) {
                        c ^= 1;
                    }
                }
            }
            q[c & 1].push_back(j);
        }

        for (int i = 0; i < 2; i++) {
            cout << "? ";
            string s(1 << n, '0');
            for (auto v : q[i]) {
                s[v] = '1';
            }
            cout << s << endl;
            ll res;
            if (!(cin >> res) || res == -1) return;
            tot += res;
        }
    }

    ll res = 2LL * ((1LL << n) - 1) - ((1 << (n - 1)) - 1);
    assert(tot % res == 0);
    cout << "! " << tot / res << endl;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
