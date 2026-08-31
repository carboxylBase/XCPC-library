#include <bits/stdc++.h>

using namespace std;
using pii = pair<int, int>;

const int N = 2000000;
const int B = sqrt(1e5);

struct Dsu {
    int fa[N];

    void init(int n) {
        for (int i = 1; i < n + 1; i++) {
            fa[i] = i;
        }
    }

    int find(int x) {
        if (x == fa[x]) {
            return x;
        }
        return fa[x] = find(fa[x]);
    }

    void merge(int x, int y) {
        x = find(x), y = find(y);
        if (x == y) {
            return;
        }
        fa[x] = y;
    }
} solver;

void solve() {
    int n, m;
    cin >> n >> m;
    vector<vector<pii>> e(m + 1);
    for (int i = 1; i <= m; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        e[c].push_back(pii(a, b));
    }

    int q;
    cin >> q;
    vector<pii> Q(q);
    map<pii, int> ans;
    for (auto &v : Q) {
        cin >> v.first >> v.second;
        ans[v] = 0;
    }

    solver.init(n + 1);

    vector<bool> vis(n + 1, false);
    vector<int> b;

    for (int i = 1; i <= m; i++) {
        for (auto [u, v] : e[i]) {
            if (!vis[u]) {
                vis[u] = true;
                b.push_back(u);
            }
            if (!vis[v]) {
                vis[v] = true;
                b.push_back(v);
            }
            solver.merge(u, v);
        }

        if (e[i].size() > B) {
            for (auto [u, v] : Q) {
                if (solver.find(u) == solver.find(v)) {
                    ans[pii(u, v)]++;
                }
            }
        } else {
            for (int x = 0; x < static_cast<int>(b.size()); x++) {
                for (int y = x + 1; y < static_cast<int>(b.size()); y++) {
                    if (solver.find(b[x]) == solver.find(b[y])) {
                        if (ans.count(pii(b[x], b[y]))) {
                            ans[pii(b[x], b[y])]++;
                        }
                        if (ans.count(pii(b[y], b[x]))) {
                            ans[pii(b[y], b[x])]++;
                        }
                    }
                }
            }
        }

        for (auto v : b) {
            solver.fa[v] = v;
            vis[v] = false;
        }
        b.clear();
    }

    for (auto v : Q) {
        cout << ans[v] << '\n';
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
    return 0;
}
