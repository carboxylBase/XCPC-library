#include <bits/stdc++.h>

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
Point<T> rotate90(const Point<T>& v) {
    return {-v.y, v.x};
}

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

// 单次 log 的复杂度
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
    int n, q; cin >> n >> q;
    vector<Point<ll>> a(n);
    map<Point<ll>, int> idx;
    for (int i = 0;i<n;i++) {
        cin >> a[i].x >> a[i].y;
        idx[a[i]] = i + 1;
    }
    vector<Point<ll>> b = getConvexHull(a);

    if (b.size() < 3) {
        while (q --) {
            int x; cin >> x >> x;
            cout << -1 << endl;
        }
        return;
    }

    Point<ll> P = b[0];
    if (b.size() == 4) {
        for (int i = 0;i<a.size();i++) {
            Point<ll> lhs = a[i] - b[0];
            Point<ll> rhs = a[i] - b[2];
            if ((lhs ^ rhs) == 0) {
                continue;
            }
            lhs = a[i] - b[1];
            rhs = a[i] - b[3];
            if ((lhs ^ rhs) != 0 ) {
                P = a[i];
                break;
            }
        }
    }

    while (q --) {
        Point<ll> z; cin >> z.x >> z.y;
        if (inConvexStrict(b, z) == 0) {
            cout << -1 << endl;
            continue;
        }
        if (b.size() == 3) {
            cout << 3 << endl;
            cout << idx[b[0]] << ' ' << idx[b[1]] << ' ' << idx[b[2]] << endl;
            continue;
        }

        bool bad = 0;
        int lo = 1, hi = b.size() - 1, nxt = 0;
        while (lo <= hi) {
            int mid = lo + hi >> 1;
            Point<ll> lhs = b[mid] - b[0];
            Point<ll> rhs = z - b[0];

            if ((lhs ^ rhs) >= 0) {
                if ((lhs ^ rhs) == 0) bad = 1;
                nxt = mid;
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }

        if (!bad) {
            cout << 3 << endl;
            cout << idx[b[0]] << ' ' << idx[b[nxt]] << ' ' << idx[b[(nxt + 1) % int(b.size())]] << endl;
            continue;
        }

        vector<Point<ll>> c;

        auto out = [&]() -> bool {
            bool bad = 0;
            for (int i = 0;i<c.size()&&!bad;i++) {
                for (int j = i + 1;j<c.size()&&!bad;j++) {
                    for (int k = j + 1;k<c.size()&&!bad;k++) {
                        vector<Point<ll>> x;
                        x.push_back(c[i]);
                        x.push_back(c[j]);
                        x.push_back(c[k]);
                        x = getConvexHull(x);
                        if (inConvexStrict(x, z)) {
                            cout << 3 << endl;
                            cout << idx[x[0]] << ' ' << idx[x[1]] << ' ' << idx[x[2]] << endl;
                            return 1;
                        }
                    }
                }
            }
            return 0;
        };

        if (b.size() > 4) {
            c.push_back(b[0]);
            c.push_back(b[nxt]);
            c.push_back(b[(nxt + 1) % int(b.size())]);
            c.push_back(b[(nxt - 1 + int(b.size())) % int(b.size())]);
            c.push_back(b[(nxt + 2) % int(b.size())]);
            c.push_back(b[(nxt - 2 + int(b.size())) % int(b.size())]);
            out();
        } else {
            c = b;
            c.push_back(P);
            if (out()) {
                continue;
            }
            cout << 4 << endl;
            for (auto v : b) {
                cout << idx[v] << ' ';
            }
            cout << endl;
        }
    }
    return;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t = 1;
    cin >> t;
    while (t--) {
        solve();
    }
    return 0;
}
