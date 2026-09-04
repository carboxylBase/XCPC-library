#include <bits/stdc++.h>

using namespace std;
using pii = pair<int, int>;

void solve() {
    int n;
    cin >> n;

    int k = 0;
    while ((1 << k) <= n) k++;
    k--;

    vector<pii> ans;

    auto build = [&](int pos) -> void {
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < (1 << k); j += (1 << (i + 1))) {
                for (int p = 0; p < (1 << i); p++) {
                    ans.emplace_back(j + p + 1 + pos,
                                     j + p + 1 + (1 << i) + pos);
                }
            }
        }
    };

    build(0);
    build(n - (1 << k));

    cout << ans.size() << '\n';
    for (auto [x, y] : ans) {
        cout << x << ' ' << y << '\n';
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
    return 0;
}
