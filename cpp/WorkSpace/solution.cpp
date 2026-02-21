#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef long long ll;


int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    if (!(cin >> n)) return 0;
    vector<Point> pts(n);
    for (int i = 0; i < n; i++) cin >> pts[i].x >> pts[i].y;

    vector<Point> hull = getConvexHull(pts);
    cout << rotatingCalipers(hull) << endl;

    return 0;
}