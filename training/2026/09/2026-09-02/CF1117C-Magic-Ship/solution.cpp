#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pll = pair<ll, ll>;

void solve() {
    ll x1, y1;
    cin >> x1 >> y1;
    ll x2, y2;
    cin >> x2 >> y2;

    int n;
    cin >> n;
    string s;
    cin >> s;

    pll period = {0, 0};
    for (char c : s) {
        if (c == 'U') {
            period.second++;
        } else if (c == 'D') {
            period.second--;
        } else if (c == 'L') {
            period.first--;
        } else {
            period.first++;
        }
    }

    ll ans = LLONG_MAX;
    pll prefix = {0, 0};

    for (int i = 0; i < n; i++) {
        ll lo = 0, hi = 100000000000LL;

        while (lo <= hi) {
            ll cycles = (lo + hi) >> 1;
            ll days = cycles * n + i;

            ll nx = x1 + period.first * cycles + prefix.first;
            ll ny = y1 + period.second * cycles + prefix.second;
            ll distance = llabs(nx - x2) + llabs(ny - y2);

            if (distance <= days) {
                ans = min(ans, days);
                hi = cycles - 1;
            } else {
                lo = cycles + 1;
            }
        }

        if (s[i] == 'U') {
            prefix.second++;
        } else if (s[i] == 'D') {
            prefix.second--;
        } else if (s[i] == 'L') {
            prefix.first--;
        } else {
            prefix.first++;
        }
    }

    cout << (ans == LLONG_MAX ? -1 : ans) << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
    return 0;
}
