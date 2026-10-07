#include <bits/stdc++.h>
using namespace std;
using ll = long long;

void solve() {
    ll x1, y1, x2, y2;
    cin >> x1 >> y1 >> x2 >> y2;
    if ((x1 == 0 && y1 == 0) || (x2 == 0 && y2 == 0)) {
        cout << "YES\n";
        return;
    }
    // Preserve the user's case analysis.
    bool ok = false;
    if (y2 == 0) {
        ok = y1 > 0 || x1 + 1 < x2;
    } else if (y1 < y2) {
        if (x1 < x2) ok = true;
        else if (x1 == x2) ok = y1 + 1 != y2;
        else ok = x2 == 0 || x1 - x2 + 1 < y2 - y1;
    } else if (y1 == y2) {
        ok = x2 == 0 || x1 + 1 < x2;
    } else {
        if (x1 + y1 - y2 + 1 < x2) ok = true;
        else if (x1 > x2) {
            ll l1 = y1 - y2 - 1 + x1 - x2 - 1;
            ll l2 = x2 - 1 + y2 - 1;
            ok = x2 == 0 || l1 >= l2;
        }
    }
    cout << (ok ? "YES\n" : "NO\n");
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    cin >> t;
    while (t--) solve();
}
