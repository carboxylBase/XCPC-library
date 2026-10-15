#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using pll = pair<ll,ll>;
void solve() {
    int n; cin >> n;
    vector<ll> a(2*n);
    for (auto &v : a) cin >> v;
    sort(a.begin(), a.end());
    vector<pll> b;
    for (int i = 0; i < 2*n; i++) {
        int j = i;
        for (; j < 2*n; j++) if (a[j] != a[i]) break;
        b.push_back(pll(a[i], j-i));
        i = j-1;
    }
    ll tot = 0;
    vector<ll> c;
    for (auto v : b) {
        if (v.second & 1) c.push_back(v.first);
        int cnt = v.second / 2;
        if (cnt & 1) tot ^= v.first;
    }
    if (c.size() > 2) { cout << "Bot\n"; return; }
    if (c.empty()) {
        if (tot == 0) cout << "Menji\n";
        else cout << "Bot\n";
    } else {
        for (auto v : c) if (v == tot) { cout << "Menji\n"; return; }
        cout << "Bot\n";
    }
}
int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int t; cin >> t;
    while (t--) solve();
}
