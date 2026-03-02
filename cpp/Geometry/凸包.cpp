#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

using db = double;
const db eps = 1e-9;

struct Point {
    double x, y;
    Point operator-(const Point& b) const { return {x - b.x, y - b.y}; }
    Point operator+(const Point& b) const { return {x + b.x, y + b.y}; }
    Point operator*(double k) const { return {x * k, y * k}; }
    double operator^(const Point& b) const { return x * b.y - y * b.x; }
    double operator&(const Point& b) const { return x * b.x + y * b.y; }
    double len2() const { return x * x + y * y; }
    double len() const { return sqrt(len2()); }
};

Point rotate90(Point &v) { return {-v.y, v.x}; }

vector<Point> getConvexHull(vector<Point>& p) {
    int n = p.size();
    sort(p.begin(), p.end(), [](Point a, Point b) {
        return a.x < b.x || (abs(a.x - b.x) < eps && a.y < b.y);
    });
    vector<Point> hull;
    for (int i = 0; i < n; ++i) {
        while (hull.size() > 1 && ((hull.back() - hull[hull.size() - 2]) ^ (p[i] - hull.back())) <= eps)
            hull.pop_back();
        hull.push_back(p[i]);
    }
    int k = hull.size();
    for (int i = n - 2; i >= 0; --i) {
        while (hull.size() > k && ((hull.back() - hull[hull.size() - 2]) ^ (p[i] - hull.back())) <= eps)
            hull.pop_back();
        hull.push_back(p[i]);
    }
    hull.pop_back();
    return hull;
}