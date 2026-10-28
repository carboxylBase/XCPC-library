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

    int c = 0;
    map<pii, int> vis;
    auto Q = [&](int i, int j) -> int {
        if (vis.count(pii(i, j))) return vis[pii(i, j)];
        c ++;
        assert(c <= 2 * n - 3);
        cout << "? " << i << ' ' << j << endl;
        int ans;
        if (!(cin >> ans) || ans == -1) exit(0);
        return vis[pii(i, j)] = ans;
    };

    int m = n;
    for (int i = 1;i<n;i++) {
        if (Q(i, i + 1) == 1) {
            m = i;
            break;
        }
    }

    if (m == n) {
        cout << "! ";
        for (int i = 1;i<n+1;i++) cout << i << ' ';
        cout << endl;
        return;
    }

    int k = m;
    for (int i = 1;i<m;i++) {
        if (Q(i, m + 1) == 1) {
            k = i;
            break;
        }
    }

    deque<int> A, B;
    for (int i = 1;i<k;i++) A.push_back(i);
    A.push_back(m + 1);
    for (int i = k;i<=m;i++) B.push_back(i);

    for (int i = m + 2;i<n+1;i++) {
        if (Q(i, A.front()) == 0) {
            A.push_front(i);
        } else if (Q(B.back(), i) == 0) {
            B.push_back(i);
        } else {
            A.push_back(i);
        }
    }

    cout << "! ";
    for (auto v : A) cout << v << ' ';
    for (auto v : B) cout << v << ' ';
    cout << endl;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t = 1;
    if (!(cin >> t)) return 0;
    while (t--) {
        solve();
    }
    return 0;
}
