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
/*
观察 getS 函数， 发现最后 tail 和 head 之间的就是半平面交的凸包
复杂度是 nlogn
*/

const db eps = 1e-10;

struct Point {
    db x, y;
    Point operator-(const Point& b) const { return {x - b.x, y - b.y}; }
    Point operator+(const Point& b) const { return {x + b.x, y + b.y}; }
    Point operator*(db k) const { return {x * k, y * k}; }
    db operator^(const Point& b) const { return x * b.y - y * b.x; }
};

struct Line {
    Point p1, p2;
    db angle;
    Line() {}
    Line(Point a, Point b) : p1(a), p2(b) {
        angle = atan2(b.y - a.y, b.x - a.x);
    }
};

Point getIntersect(Line a, Line b) {
    Point v1 = a.p2 - a.p1;
    Point v2 = b.p2 - b.p1;
    Point u = a.p1 - b.p1;
    db t = (v2 ^ u) / (v1 ^ v2);
    return a.p1 + v1 * t;
}

bool isRight(Line l, Point p) {
    return ((l.p2 - l.p1) ^ (p - l.p1)) < -eps;
}

db getS(vector<Line>& lines) {
    sort(lines.begin(), lines.end(), [](Line a, Line b) {
        if (abs(a.angle - b.angle) > eps) return a.angle < b.angle;
        return ((a.p2 - a.p1) ^ (b.p2 - a.p1)) < -eps;
    });

    int n = lines.size();
    int m = 0;
    for (int i = 0; i < n; i++) {
        if (i > 0 && abs(lines[i].angle - lines[i - 1].angle) < eps) continue;
        lines[m++] = lines[i];
    }
    lines.erase(lines.begin() + m, lines.end());

    int head = 0, tail = 0;
    vector<Line> dq(lines.size() + 5);
    vector<Point> p(lines.size() + 5);

    dq[tail++] = lines[0];
    for (int i = 1; i < (int)lines.size(); i++) {
        while (tail - head > 1 && isRight(lines[i], p[tail - 1])) tail--;
        while (tail - head > 1 && isRight(lines[i], p[head + 1])) head++;
        
        dq[tail++] = lines[i];
        if (abs((dq[tail - 1].p2 - dq[tail - 1].p1) ^ (dq[tail - 2].p2 - dq[tail - 2].p1)) < eps) {
            tail--;
            if (!isRight(dq[tail], lines[i].p1)) dq[tail] = lines[i];
        }
        if (tail - head > 1) p[tail - 1] = getIntersect(dq[tail - 2], dq[tail - 1]);
    }

    while (tail - head > 1 && isRight(dq[head], p[tail - 1])) tail--;
    while (tail - head > 1 && isRight(dq[tail - 1], p[head + 1])) head++;

    if (tail - head < 3) return 0.0;

    p[head] = getIntersect(dq[tail - 1], dq[head]);

    db area = 0;
    for (int i = head; i < tail; i++) {
        area += (p[i] ^ p[i == tail - 1 ? head : i + 1]);
    }
    return abs(area) / 2.0;
}