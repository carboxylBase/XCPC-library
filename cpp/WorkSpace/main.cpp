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

namespace Rotating_Calipers {

struct Point {
    ll x, y;
    Point operator-(const Point& b) const { return {x - b.x, y - b.y}; }
    ll operator^(const Point& b) const { return x * b.y - y * b.x; }
};

ll distSq(Point a, Point b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

ll cross_product(Point a, Point b, Point c) {
    return (b - a) ^ (c - a);
}

vector<Point> getConvexHull(vector<Point>& pts) {
    int n = pts.size();
    if (n <= 2) return pts;
    sort(pts.begin(), pts.end(), [](Point a, Point b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    vector<Point> hull;
    for (int i = 0; i < n; ++i) {
        while (hull.size() >= 2 && cross_product(hull[hull.size() - 2], hull.back(), pts[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(pts[i]);
    }
    int lower_size = hull.size();
    for (int i = n - 2; i >= 0; --i) {
        while (hull.size() > lower_size && cross_product(hull[hull.size() - 2], hull.back(), pts[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(pts[i]);
    }
    hull.pop_back();
    return hull;
}

ll rotatingCalipers(const vector<Point>& hull) {
    int m = hull.size();
    if (m == 2) return distSq(hull[0], hull[1]);

    ll max_dist = 0;
    int j = 1;

    for (int i = 0; i < m; i++) {
        int next_i = (i + 1) % m;
        while (abs(cross_product(hull[i], hull[next_i], hull[(j + 1) % m])) > 
               abs(cross_product(hull[i], hull[next_i], hull[j]))) {
            j = (j + 1) % m;
        }
        max_dist = max({max_dist, distSq(hull[i], hull[j]), distSq(hull[next_i], hull[j])});
    }
    return max_dist;
}

};
using namespace Rotating_Calipers;

void solve() {
	int n; cin >> n;
	vector<Point> a(n);
	for (int i = 0;i<n;i++) cin >> a[i].x >> a[i].y;
	vector<Point> b = getConvexHull(a);
	cout << rotatingCalipers(b) << endl;
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
