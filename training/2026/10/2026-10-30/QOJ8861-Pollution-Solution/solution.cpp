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
const db PI = acosl(-1);

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

using P = Point<db>;

vector<P> splitEdge(P A, P B, P C, db r,
                    db start, db sweep) {
    const db TAU = 2 * acosl(-1.0L);
    const db ANG_EPS = 1e-12L;

    assert(r > 0);
    assert(0 <= sweep && sweep <= TAU);

    // 归一化到 [0, 2π)
    auto norm = [&](db angle) {
        angle = fmodl(angle, TAU);
        if (angle < 0) angle += TAU;
        return angle;
    };

    start = norm(start);

    // v 是相对于圆心的向量
    auto onArc = [&](const P &v) {
        db rel = norm(atan2l(v.y, v.x) - start);

        // 第二项处理起始角附近的浮点误差
        return rel <= sweep + ANG_EPS ||
               rel >= TAU - ANG_EPS;
    };

    P u = A - C, d = B - A;
    db a = d & d;
    if (a == 0) return {A};

    // 圆心到直线 AB 的垂足，在线段参数中的位置
    db t0 = -(u & d) / a;

    db cr = d ^ u;
    db dist2 = cr * cr / a;
    db h2 = r * r - dist2;

    db tol = 32 * numeric_limits<db>::epsilon()
             * max({1.0L, r * r, dist2});

    if (h2 < -tol) return {A, B};

    h2 = max(0.0L, h2);
    db dt = sqrtl(h2 / a);

    // 垂足相对于圆心的坐标
    P foot = P{-d.y, d.x} * (cr / a);

    vector<pair<db, P>> cuts{{0, A}, {1, B}};

    auto add = [&](db t, P v) {
        if (0 < t && t < 1 && onArc(v)) {
            cuts.push_back({t, C + v});
        }
    };

    add(t0 - dt, foot - d * dt);
    add(t0 + dt, foot + d * dt);

    sort(cuts.begin(), cuts.end(),
         [](const auto &x, const auto &y) {
             return x.first < y.first;
         });

    cuts.erase(unique(cuts.begin(), cuts.end(),
                      [](const auto &x, const auto &y) {
                          return x.first == y.first;
                      }),
               cuts.end());

    vector<P> res;
    for (auto &[t, p] : cuts) res.push_back(p);
    return res;
}

void solve() {
    int n, r; cin >> n >> r;
    vector<P> a(n), q;
    for (auto &v : a) cin >> v.x >> v.y;
    for (int i = 0;i<n;i++) {
        P Z;
        Z.x = Z.y = 0;
        auto seg = splitEdge(a[i], a[(i + 1) % n], Z, r, 0, PI);
        for (int j = 0;j+1<seg.size();j++) {
            q.push_back(seg[j]);
        }
    }

    db ans = 0;
    for (int i = 0;i<q.size();i++) {
        int j = (i + 1) % int(q.size());
        db x = (q[i].x + q[j].x) / 2;
        db y = (q[i].y + q[j].y) / 2;
        if (x * x + y * y > r * r) {
            db val_c = (q[i] ^ q[j]);
            db val_d = (q[i] & q[j]);
            db cur = atan2l(val_c, val_d);
            ans += r * r * cur / 2;
        } else {
            ans += (q[i] ^ q[j]) / 2;
        }
    }
    cout << fixed << setprecision(12) << ans << endl;
    return;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t = 1;
    // cin >> t;

    while (t--) {
        solve();
    }

    return 0;
}
