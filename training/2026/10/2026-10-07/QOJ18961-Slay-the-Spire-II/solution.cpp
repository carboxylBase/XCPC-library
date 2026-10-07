#include <bits/stdc++.h>
using namespace std;
using ll = long long;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<ll> a, b;
    ll mx = 0;
    for (int i = 0; i < n; i++) {
        char op;
        ll d;
        cin >> op >> d;
        if (op == '1') a.push_back(d);
        else if (op == '2') b.push_back(d);
        else mx = max(mx, d);
    }
    sort(a.begin(), a.end(), greater<ll>());
    sort(b.begin(), b.end(), greater<ll>());
    int p = -1;
    for (int i = 0; i < (int)b.size(); i++) {
        if (b[i] >= 2 * mx) p = i;
    }
    for (int i = 1; i < (int)a.size(); i++) a[i] += a[i - 1];
    for (int i = 1; i < (int)b.size(); i++) b[i] += b[i - 1];
    ll ans = mx * m;
    for (int i = 0; i <= min(m, (int)a.size()); i++) {
        int pos = min((m - i) / 2, p + 1) - 1;
        ll res = 0;
        if (i > 0) res += a[i - 1];
        if (pos >= 0) res += b[pos];
        res += (m - ((pos + 1) * 2 + i)) * mx;
        ans = max(ans, res);
    }
    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
