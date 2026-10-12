#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using db = double;
using i128 = __int128_t;

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

template<class T>
vector<Point<T>> minkowski(vector<Point<T>> a, vector<Point<T>> b) {
    if (a.empty() || b.empty()) return {};

    if (a.size() == 1) {
        for (auto& v : b) v = v + a[0];
        return b;
    }
    if (b.size() == 1) {
        for (auto& v : a) v = v + b[0];
        return a;
    }

    auto normalize = [](vector<Point<T>>& v) {
        auto it = min_element(v.begin(), v.end(),
            [](const Point<T>& a, const Point<T>& b) {
                return a.y < b.y || (a.y == b.y && a.x < b.x);
            });
        rotate(v.begin(), it, v.end());
    };

    auto sgn = [](T v) -> int {
        if constexpr (is_floating_point_v<T>) {
            return (v > eps) - (v < -eps);
        } else {
            return (v > 0) - (v < 0);
        }
    };

    normalize(a);
    normalize(b);

    int n = a.size(), m = b.size();
    vector<Point<T>> va(n), vb(m);

    for (int i = 0; i < n; i++) {
        va[i] = a[(i + 1) % n] - a[i];
    }
    for (int i = 0; i < m; i++) {
        vb[i] = b[(i + 1) % m] - b[i];
    }

    vector<Point<T>> ans;
    ans.reserve(n + m + 1);
    ans.push_back(a[0] + b[0]);

    int i = 0, j = 0;
    while (i < n && j < m) {
        int s = sgn(va[i] ^ vb[j]);
        if (s == 0) {
            ans.push_back(ans.back() + va[i++] + vb[j++]);
        } else if (s > 0) {
            ans.push_back(ans.back() + va[i++]);
        } else {
            ans.push_back(ans.back() + vb[j++]);
        }
    }

    while (i < n) ans.push_back(ans.back() + va[i++]);
    while (j < m) ans.push_back(ans.back() + vb[j++]);

    ans.pop_back();
    return ans;
}