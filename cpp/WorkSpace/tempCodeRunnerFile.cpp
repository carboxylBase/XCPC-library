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

ll dp[101][2001];
ll nxt[101][2001];

void solve() {
	int n, x, y, t; cin >> n >> x >> y >> t;
	vector<ll> a(n + 1, 0);
	for (int i = 1;i<n+1;i++) cin >> a[i];
	vector<ll> b(n + 1, 0);
	for (int i = 1;i<n+1;i++) cin >> b[i];
	for (int i = 1;i<=n;i++) {
		for (int j = 0;j<=t;j++) {
			for (int k = 0;k<=a[i];k++) {
				if (k && x + y * k + j <= t) {
					if (dp[i][x + y * k + j] < dp[i - 1][j] + k * b[i]) {
						nxt[i][x + y * k + j] = k;
					}
					dp[i][x + y * k + j] = max(dp[i][x + y * k + j], dp[i - 1][j] + k * b[i]);
				}
				if (!k && j <= t) {
					if (dp[i][j] < dp[i - 1][j]) {
						nxt[i][j] = 0;
					}
					dp[i][j] = max(dp[i][j], dp[i - 1][j]);
				}
			}
		}
	}

	ll mx = 0;
	int tot = 0;
	for (int i = 0;i<=t;i++) {
		if (mx < dp[n][i]) tot = i;
		mx = max(mx, dp[n][i]);
	}

	vector<int> ans(n + 1, 0);
	for (int i = n;i>=1;i--) {
		ans[i] = tot;
		tot -= nxt[i][tot] * b[i];
	}

	for (int i = 1;i<n+1;i++) {
		cout << ans[i] << ' ';
	}
	cout << endl;
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
