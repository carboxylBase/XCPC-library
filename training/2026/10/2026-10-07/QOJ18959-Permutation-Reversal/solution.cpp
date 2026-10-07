#include <bits/stdc++.h>
using namespace std;
using ll = long long;

void solve() {
    int n;
    cin >> n;
    vector<ll> a(n + 1, 0);
    for (int i = 1; i <= n; i++) cin >> a[i];
    vector<int> s(1, 0), ls(n + 1, 0), rs(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        while (a[s.back()] > a[i]) {
            ls[i] = s.back();
            s.pop_back();
        }
        rs[s.back()] = i;
        s.push_back(i);
    }
    int root = rs[0];
    // Iterative traversal avoids recursion depth n on a monotone permutation.
    vector<int> order{root}, left(n + 1), right(n + 1);
    left[root] = 1;
    right[root] = n;
    for (int i = 0; i < (int)order.size(); i++) {
        int u = order[i];
        if (ls[u]) {
            left[ls[u]] = left[u];
            right[ls[u]] = u - 1;
            order.push_back(ls[u]);
        }
        if (rs[u]) {
            left[rs[u]] = u + 1;
            right[rs[u]] = right[u];
            order.push_back(rs[u]);
        }
    }
    vector<ll> mn(n + 1), mx(n + 1), sum(n + 1);
    for (int i = n - 1; i >= 0; i--) {
        int u = order[i], l = left[u], r = right[u];
        sum[u] = a[u] + sum[ls[u]] + sum[rs[u]];
        mn[u] = mn[ls[u]] + mn[rs[u]] + u * a[u];
        mx[u] = mx[ls[u]] + mx[rs[u]] + u * a[u];
        ll shifted = -(u + 1 - l) * sum[rs[u]]
                     + (r - u + 1) * sum[ls[u]]
                     + (l + r - u) * a[u];
        mn[u] = min(mn[u], mn[ls[u]] + mn[rs[u]] + shifted);
        mx[u] = max(mx[u], mx[ls[u]] + mx[rs[u]] + shifted);
    }
    cout << mx[root] << ' ' << mn[root] << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
