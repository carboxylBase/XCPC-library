// 根据用户实现与用户确认 AC 的卡壳修正整理，非原始提交文件的逐字副本。
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using db = long double;
const ll INF = 5e18;
const db eps = 1e-9;

template<class T>
struct Point {
    T x, y;

    Point operator-(const Point& b) const {
        return {x - b.x, y - b.y};
    }

    Point operator+(const Point& b) const {
        return {x + b.x, y + b.y};
    }

    Point operator*(T k) const {
        return {x * k, y * k};
    }

    T operator^(const Point& b) const {
        return x * b.y - y * b.x;
    }

    T operator&(const Point& b) const {
        return x * b.x + y * b.y;
    }

    bool operator<(const Point& b) const {
        return x < b.x || (x == b.x && y < b.y);
    }

    bool operator==(const Point& b) const {
        return x == b.x && y == b.y;
    }

    T len2() const {
        return x * x + y * y;
    }

    db len() const {
        return hypot((db)x, (db)y);
    }
};

template<class T>
vector<Point<T>> getConvexHull(vector<Point<T>>& p) {
    sort(p.begin(), p.end());
    p.erase(unique(p.begin(), p.end()), p.end());

    int n = p.size();
    if (n <= 1) return p;

    auto check = [](T v) {
        if constexpr (is_floating_point_v<T>) {
            return v <= eps;
        } else {
            return v <= 0;
        }
    };

    vector<Point<T>> hull;

    for (int i = 0; i < n; i++) {
        while (hull.size() > 1 &&
               check((hull.back() - hull[hull.size() - 2]) ^
                     (p[i] - hull.back()))) {
            hull.pop_back();
        }
        hull.push_back(p[i]);
    }

    int k = hull.size();
    for (int i = n - 2; i >= 0; i--) {
        while ((int)hull.size() > k &&
               check((hull.back() - hull[hull.size() - 2]) ^
                     (p[i] - hull.back()))) {
            hull.pop_back();
        }
        hull.push_back(p[i]);
    }

    hull.pop_back();
    return hull;
}

template<class T>
bool inConvexStrict(const vector<Point<T>>& hull, const Point<T>& p) {
    int n = hull.size();
    if (n < 3) return false;

    auto cross = [](const Point<T>& a,
                    const Point<T>& b,
                    const Point<T>& c) -> T {
        return (b - a) ^ (c - a);
    };

    if (cross(hull[0], hull[1], p) <= 0) return false;
    if (cross(hull[0], hull[n - 1], p) >= 0) return false;

    int lo = 1, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (cross(hull[0], hull[mid], p) >= 0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return cross(hull[lo], hull[hi], p) > 0;
}

void solve() {
    int n; cin >> n;
    vector<Point<ll>> a(n);
    for (auto &v : a) cin >> v.x >> v.y;
    vector<Point<ll>> b = getConvexHull(a);

    if (a.size() == b.size()) {
        cout << -1 << endl;
        return;
    }

    vector<Point<ll>> c;
    for (auto v : a) {
        if (inConvexStrict(b, v)) c.push_back(v);
    }
    c = getConvexHull(c);

    int h = b.size(), m = c.size();
    auto area = [&](int i, int j) -> ll {
        return (b[(i + 1) % h] - b[i]) ^ (c[j] - b[i]);
    };

    int j = 0;
    for (int t = 1; t < m; ++t) {
        if (area(0, t) < area(0, j)) j = t;
    }

    ll ans = INF;
    for (int i = 0; i < h; ++i) {
        while (area(i, (j + 1) % m) < area(i, j)) {
            j = (j + 1) % m;
        }
        ans = min(ans, area(i, j));
    }

    ans = -ans;
    ll res = 0;
    for (int i = 2; i < (int)b.size(); i++) {
        res += abs((b[i] - b[0]) ^ (b[i - 1] - b[0]));
    }
    ans += abs(res);
    cout << ans << endl;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t = 1;
    cin >> t;
    while (t--) solve();
    return 0;
}
