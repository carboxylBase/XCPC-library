#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pll = pair<ll, ll>;

void solve() {
    int n;
    cin >> n;

    vector<pll> a(n);
    for (auto &v : a) cin >> v.first >> v.second;

    for (int i = n - 1; i >= 0; i--) {
        a[i].first -= a[0].first;
        a[i].second -= a[0].second;
    }

    // 00, 01, 10, 11: the parities of x and y.
    while (true) {
        int cnt[4] = {0, 0, 0, 0};
        for (auto [x, y] : a) {
            int type;
            if (x % 2 == 0) {
                type = (y % 2 == 0 ? 0 : 1);
            } else {
                type = (y % 2 == 0 ? 2 : 3);
            }
            cnt[type]++;
        }

        if (cnt[0] == n || cnt[1] == n ||
            cnt[2] == n || cnt[3] == n) {
            for (auto &v : a) {
                v.first /= 2;
                v.second /= 2;
            }
        } else {
            break;
        }
    }

    vector<int> q[4];
    for (int i = 0; i < n; i++) {
        auto [x, y] = a[i];
        int type;
        if (x % 2 == 0) {
            type = (y % 2 == 0 ? 0 : 1);
        } else {
            type = (y % 2 == 0 ? 2 : 3);
        }
        q[type].push_back(i);
    }

    vector<int> ans;
    if (q[0].size() + q[3].size() != static_cast<size_t>(n)) {
        ans.insert(ans.end(), q[0].begin(), q[0].end());
        ans.insert(ans.end(), q[3].begin(), q[3].end());
    } else {
        ans = q[0];
    }

    cout << ans.size() << '\n';
    for (int v : ans) cout << v + 1 << ' ';
    cout << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
    return 0;
}
