#include <bits/stdc++.h>

using namespace std;

struct CustomHash {
    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    size_t operator()(uint64_t x) const {
        static const uint64_t seed =
            chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + seed);
    }
};

struct DSU {
    vector<int> parent;
    vector<int> size;

    explicit DSU(int n) : parent(n + 1), size(n + 1) {}

    void init(const vector<int>& vertices) {
        for (int v : vertices) {
            parent[v] = v;
            size[v] = 1;
        }
    }

    int find(int x) {
        if (parent[x] == x) return x;
        return parent[x] = find(parent[x]);
    }

    void merge(int x, int y) {
        x = find(x);
        y = find(y);
        if (x == y) return;

        if (size[x] < size[y]) swap(x, y);
        parent[y] = x;
        size[x] += size[y];
    }
};

struct Edge {
    int u;
    int v;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n, m;
    cin >> n >> m;

    vector<vector<Edge>> edgesByColor(m + 1);
    for (int i = 0; i < m; ++i) {
        int u, v, color;
        cin >> u >> v >> color;
        edgesByColor[color].push_back({u, v});
    }

    int q;
    cin >> q;

    auto makeKey = [n](int u, int v) -> long long {
        if (u > v) swap(u, v);
        return 1LL * u * (n + 1) + v;
    };

    unordered_map<long long, int, CustomHash> queryIndex;
    queryIndex.reserve(2 * q + 1);
    queryIndex.max_load_factor(0.7f);

    vector<pair<int, int>> uniqueQueries;
    vector<int> queryId(q);

    for (int i = 0; i < q; ++i) {
        int u, v;
        cin >> u >> v;
        if (u > v) swap(u, v);

        long long key = makeKey(u, v);
        auto [it, inserted] = queryIndex.emplace(
            key, static_cast<int>(uniqueQueries.size())
        );

        if (inserted) uniqueQueries.emplace_back(u, v);
        queryId[i] = it->second;
    }

    vector<int> answer(uniqueQueries.size(), 0);
    vector<int> mark(n + 1, 0);
    DSU dsu(n);

    const int blockSize = static_cast<int>(sqrt(m)) + 1;
    int stamp = 0;

    for (int color = 1; color <= m; ++color) {
        const auto& edges = edgesByColor[color];
        if (edges.empty()) continue;

        vector<int> vertices;
        vertices.reserve(2 * edges.size());
        for (const auto& [u, v] : edges) {
            vertices.push_back(u);
            vertices.push_back(v);
        }

        sort(vertices.begin(), vertices.end());
        vertices.erase(unique(vertices.begin(), vertices.end()), vertices.end());

        ++stamp;
        dsu.init(vertices);
        for (int v : vertices) mark[v] = stamp;
        for (const auto& [u, v] : edges) dsu.merge(u, v);
        for (int v : vertices) dsu.parent[v] = dsu.find(v);

        if (static_cast<int>(edges.size()) <= blockSize) {
            vector<pair<int, int>> componentVertices;
            componentVertices.reserve(vertices.size());

            for (int v : vertices) {
                componentVertices.emplace_back(dsu.parent[v], v);
            }
            sort(componentVertices.begin(), componentVertices.end());

            for (int left = 0; left < static_cast<int>(componentVertices.size());) {
                int right = left;
                while (right < static_cast<int>(componentVertices.size()) &&
                       componentVertices[right].first == componentVertices[left].first) {
                    ++right;
                }

                for (int i = left; i < right; ++i) {
                    for (int j = i + 1; j < right; ++j) {
                        int u = componentVertices[i].second;
                        int v = componentVertices[j].second;
                        long long key = makeKey(u, v);

                        auto it = queryIndex.find(key);
                        if (it != queryIndex.end()) ++answer[it->second];
                    }
                }

                left = right;
            }
        } else {
            for (int i = 0; i < static_cast<int>(uniqueQueries.size()); ++i) {
                auto [u, v] = uniqueQueries[i];
                if (mark[u] == stamp && mark[v] == stamp &&
                    dsu.parent[u] == dsu.parent[v]) {
                    ++answer[i];
                }
            }
        }
    }

    for (int id : queryId) {
        cout << answer[id] << '\n';
    }

    return 0;
}
