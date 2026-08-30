#define DEBUG 1
#define FUCK cout << "fuck" << endl;
#if DEBUG
	#include "all.hpp"
#else
	#include <bits/stdc++.h>
#endif

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
	ll n, m, k; cin >> n >> m >> k;
	vector<ll> a(n);
	for (auto &v : a) cin >> v;
	while (m --) {
		char opt; cin >> opt;
		if (opt != 'A') {
			ll t; cin >> t;
			while (t --) {
				ll mx = -INF;
				for (auto v : a) mx = max(v, mx);
				for (auto &v : a) {
					if (v == mx) {
						v -= k;
						break;
					}
				}
			}

			for (auto v : a) {
				cout << v << ' ';
			}
			cout << endl;
		} else {
		}
	}
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
	// cin >> t;

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
