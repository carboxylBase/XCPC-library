#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

using namespace std;

struct DSU {
    vector<int> parent;
    vector<int> size;

    explicit DSU(int n) : parent(n + 1), size(n + 1, 1) {
        iota(parent.begin(), parent.end(), 0);
    }

    int find(int x) {
        if (parent[x] == x) return x;
        return parent[x] = find(parent[x]);
    }

    bool merge(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) return false;
        if (size[x] < size[y]) swap(x, y);
        parent[y] = x;
        size[x] += size[y];
        return true;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m1, m2;
    cin >> n >> m1 >> m2;

    DSU dsu1(n), dsu2(n);
    for (int i = 0; i < m1; ++i) {
        int u, v;
        cin >> u >> v;
        dsu1.merge(u, v);
    }
    for (int i = 0; i < m2; ++i) {
        int u, v;
        cin >> u >> v;
        dsu2.merge(u, v);
    }

    vector<pair<int, int>> answer;

    for (int i = 2; i <= n; ++i) {
        if (dsu1.find(i) != dsu1.find(1) &&
            dsu2.find(i) != dsu2.find(1)) {
            answer.emplace_back(1, i);
            dsu1.merge(1, i);
            dsu2.merge(1, i);
        }
    }

    vector<int> left, right;
    for (int i = 2; i <= n; ++i) {
        if (dsu2.find(i) == i && dsu2.find(i) != dsu2.find(1)) {
            left.push_back(i);
        }
        if (dsu1.find(i) == i && dsu1.find(i) != dsu1.find(1)) {
            right.push_back(i);
        }
    }

    int pairs = static_cast<int>(min(left.size(), right.size()));
    for (int i = 0; i < pairs; ++i) {
        answer.emplace_back(left[i], right[i]);
        dsu1.merge(left[i], right[i]);
        dsu2.merge(left[i], right[i]);
    }

    cout << answer.size() << '\n';
    for (auto [u, v] : answer) {
        cout << u << ' ' << v << '\n';
    }
}
