#include <bits/stdc++.h>
#define DEBUG 1
using ll = long long;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;
using namespace std;

void solve() {
    int n,m; cin >> n >> m;
    vector<ll> a(n);
    for (int i = 0;i<n;i++) {
        cin >> a[i];
    }
    sort(a.begin(),a.end());

    auto check = [&](int x) -> bool {
        multiset<ll> q;
        for (int i = 0;i<x;i++) {
            q.insert(m);
        }

        for (int i = 0;i<n;i++) {
            if (*q.rbegin() < a[i]) {
                return 0;
            }
            ll z = *q.rbegin();
            q.erase(q.find(z));
            q.insert(z - a[i]);
        }

        return 1;
    };

    int ls = 1,rs = n,mid,ans;
    while (ls <= rs) {
        mid = ls + rs >> 1;
        if (check(mid)) {
            rs = mid - 1;
            ans = mid;
        } else {
            ls = mid + 1;
        }
    }

    cout << ans << endl;
    return;
}

signed main() {
    #if DEBUG
        freopen("input.txt", "r", stdin);
        auto start_time = chrono::steady_clock::now();
    #else
        ios::sync_with_stdio(false);
    #endif
    cin.tie(nullptr);

    int t = 1;
    cin >> t;

    while (t--) {
        solve();
    }

    #if DEBUG
        auto end_time = chrono::steady_clock::now();
        auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
        cerr << "Time: " << diff.count() << " ms" << endl;

    #endif

    return 0;
}