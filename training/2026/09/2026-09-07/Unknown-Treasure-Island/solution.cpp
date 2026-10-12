#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using pll = pair<ll, ll>;
using pii = pair<int, int>;
using db = long double;
using pdd = pair<db, db>;

const int N = 2e6;
const ll INF = 1e18;
const ll MOD = 1e9 + 7;

void solve() {
    ll n, k; cin >> n >> k;
    multiset<pll> a, b;
    ll A = 0, B = 0;
    for (int i = 0; i < n; i++) {
        ll o, x; cin >> o >> x;
        if (o == 1) {
            A += x;
            a.insert(pll(i, x));
        } else if (o == 2) {
            B += x;
        } else {
            ll res = min(A, x / k);
            A -= res;
            x -= res * k;

            while (res) {
                if (a.rbegin()->second <= res) {
                    res -= a.rbegin()->second;
                    a.erase(a.find(*a.rbegin()));
                } else {
                    auto it = *a.rbegin();
                    a.erase(a.find(it));
                    it.second -= res;
                    if (it.second > 0) {
                        a.insert(it);
                    }
                    res = 0;
                }
            }

            if (x <= B) {
                B -= x;
                b.insert(pll(k - x, i));
                continue;
            }

            if (A == 0) {
                for (int j = i + 1; j < n; j++) cin >> o >> x;
                cout << -1 << endl;
                return;
            }

            bool ok = 0;
            while (1) {
                if (b.empty()) {
                    break;
                }

                int pos = b.begin()->second;
                ll val = b.begin()->first;

                if (val > k - x) {
                    break;
                } else {
                    b.erase(b.begin());
                    auto it = a.lower_bound(pll(pos, -1));
                    if (!a.empty() && it != a.begin()) {
                        --it;
                        pll z = *it;
                        B += k - val;
                        B -= x;
                        A--;
                        b.insert(pll(k - x, i));
                        a.erase(it);
                        z.second--;
                        if (z.second > 0) {
                            a.insert(z);
                        }
                        ok = 1;
                        break;
                    }
                }
            }

            if (!ok) {
                A--;
                auto it = *a.rbegin();
                a.erase(a.find(it));
                it.second--;
                if (it.second > 0) {
                    a.insert(it);
                }
            }
        }
    }

    cout << A * k + B << endl;
}

signed main() {
    ios::sync_with_stdio(0); cin.tie(0);

    int _ = 1;
    cin >> _;
    while (_--) {
        solve();
    }
    return 0;
}
