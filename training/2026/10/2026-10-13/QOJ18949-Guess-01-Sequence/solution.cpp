#include <bits/stdc++.h>
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

void solve() {
    int n; cin >> n;
    int z = ceil(log2(n));

    auto Q = [&](vector<int>& q) -> ll {
        cout << "? ";
        cout << q.size() << ' ';
        for (auto v : q) cout << v << ' ';
        cout << endl;
        ll x; cin >> x;
        return x;
    };

    for (int i = z - 1;i>=0;i--) {
        vector<int> q;
        for (int j = 0;j<n;j+=(1<<(i+1))) {
            for (int k = 0;k<(1<<i);k++) {
                if (j + k + 1 <= n) {
                    q.push_back(j + k + 1);
                }
            }
        }

        ll x = Q(q);
        if (x != 0) {
            int lo = 1, hi = (1 << i) - 1;
            ll ans = x;
            while (lo <= hi) {
                int mid = lo + hi >> 1;
                vector<int> q;
                for (int j = 0;j<n;j+=(1<<(i+1))) {
                    for (int k = 0;k<mid;k++) {
                        if (j + k + 1 <= n) {
                            q.push_back(j + k + 1);
                        }
                    }
                }
                ll x = Q(q);
                if (x) {
                    ans = x;
                    hi = mid - 1;
                } else {
                    lo = mid + 1;
                }
            }
            cout << "! " << ans + 1 << endl;
            return;
        }
    }

    cout << "! 1" << endl;
    return;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t = 1;
    while (t--) solve();
    return 0;
}
