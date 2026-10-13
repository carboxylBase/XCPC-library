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

int fa[N];
ll val[N];

int F(int x) {
    if (x == fa[x]) return x;
    return fa[x] = F(fa[x]);
}

void merge(int x, int y) {
    x = F(x), y = F(y);
    fa[x] = y;
    return;
}

struct Node {
    int mn, mx, sum;
    Node() {
        mn = mx = sum = 0;
    }
} nodes[N * 4];

Node merge(Node x, Node y) {
    Node z;
    z.mn = INT_MAX, z.mx = INT_MIN;
    if (x.mn) z.mn = min(z.mn, x.mn);
    if (y.mn) z.mn = min(z.mn, y.mn);
    if (x.mx) z.mx = max(z.mx, x.mx);
    if (y.mx) z.mx = max(z.mx, y.mx);
    if (z.mn == INT_MAX) z.mn = 0;
    if (z.mx == INT_MIN) z.mx = 0;
    z.sum = x.sum + y.sum;
    return z; 
}

#define LS (rt << 1)
#define RS (rt << 1 | 1)

void modify(int rt, int l, int r, int q, int val) {
    if (l == r) {
        nodes[rt].mn = nodes[rt].mx = val;
        if (val == 0) nodes[rt].sum = 0;
        else nodes[rt].sum = 1;
        return;
    }
    int mid = l + r >> 1;
    if (q <= mid) {
        modify(LS, l, mid, q, val);
    } else {
        modify(RS, mid + 1, r, q, val);
    }
    nodes[rt] = merge(nodes[LS], nodes[RS]);
    return;
}

pii query(int rt, int l, int r, int& k) {
    if (!k) return pii(0, 0);
    if (nodes[rt].sum <= k) {
        k -= nodes[rt].sum;
        return pii(nodes[rt].mn, nodes[rt].mx);
    }
    int mid = l + r >> 1;
    pii ans = pii(INT_MAX, INT_MIN);
    pii res = query(LS, l, mid, k);
    if (res.first) ans.first = min(ans.first, res.first);
    if (res.second) ans.second = max(ans.second, res.second);
    res = query(RS, mid + 1, r, k);
    if (res.first) ans.first = min(ans.first, res.first);
    if (res.second) ans.second = max(ans.second, res.second);
    return ans;
}

void solve() {
    ll n, m, q; cin >> n >> m >> q;
    vector<int> a(n + 1, 0);
    for (int i = 1;i<n+1;i++) {
        cin >> a[i];
    }

    vector<array<ll, 3>> g(m);
    for (auto &v : g) cin >> v[1] >> v[2] >> v[0];
    sort(g.begin(), g.end());

    for (int i = 1;i<n+1;i++) {
        fa[i] = i;
    }
    int tp = n + 1;

    vector<vector<int>> e(n + 1);

    for (auto [w, u, v] : g) {
        if (F(u) == F(v)) continue;
        u = F(u), v = F(v);
        fa[tp] = tp;
        e.push_back(vector<int>());
        e[tp].push_back(u), e[tp].push_back(v);
        merge(u, tp), merge(v, tp);
        val[tp] = w;
        tp ++;
    }

    vector<pii> b;
    for (int i = 1;i<=n;i++) {
        b.push_back(pii(a[i], i));
    }

    vector<array<int, 3>> Q(q);
    for (auto &v : Q) {
        cin >> v[0];
        if (v[0] == 1) {
            cin >> v[1] >> v[2];
            b.push_back(pii(v[2], v[1]));
        } else {
            cin >> v[1];
        }
    }

    sort(b.begin(), b.end());
    b.erase(unique(b.begin(), b.end()), b.end());
    sort(b.begin(), b.end(), [&](pii l, pii r) {
        if (l.first != r.first) return l.first > r.first;
        return l.second < r.second;
    });

    auto idx = [&](pii z) -> int {
        int lo = 0, hi = b.size() - 1;
        while (lo <= hi) {
            int mid = lo + hi >> 1;
            if (z.first == b[mid].first) {
                if (z.second == b[mid].second) {
                    return mid;
                } else if (z.second > b[mid].second) {
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            } else if (z.first < b[mid].first) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
    };

    vector<int> dfn(tp, 0), dp(tp + 1, 0);
    int cnt_dfn = 0;
    int M = 22;
    vector<vector<int>> fa(M, vector<int>(tp, 0));
    vector<int> dep(tp, 0);
    dep[tp - 1] = 1;
    auto dfs = [&](auto&&self, int u) -> void {
        dfn[u] = ++ cnt_dfn;
        dp[cnt_dfn] = u;
        for (auto v : e[u]) {
            dep[v] = dep[u] + 1;
            fa[0][v] = u;
            self(self, v);
        }
        return;
    };
    dfs(dfs, tp - 1);

    for (int i = 1;i<M;i++) {
        for (int j = 1;j<tp;j++) {
            fa[i][j] = fa[i - 1][fa[i - 1][j]];
        }
    }

    auto lca = [&](int x, int y) -> int {
        for (int i = M - 1;i>=0;i--) {
            if (dep[x] > dep[y]) swap(x, y);
            if (dep[fa[i][y]] >= dep[x]) {
                y = fa[i][y];
            }
        }
        if (x == y) return x;
        for (int i = M - 1;i>=0;i--) {
            if (fa[i][x] != fa[i][y]) {
                x = fa[i][x], y = fa[i][y];
            }
        }
        return fa[0][x];
    };

    for (int i = 1;i<n+1;i++) {
        int id = idx(pii(a[i], i));
        modify(1, 0, b.size() - 1, id, dfn[i]);
    }

    for (auto v : Q) {
        if (v[0] == 1) {
            int c = v[1], t = v[2];
            int id = idx(pii(a[c], c));
            modify(1, 0, b.size() - 1, id, 0);
            a[c] = t;
            id = idx(pii(a[c], c));
            modify(1, 0, b.size() - 1, id, dfn[c]);
        } else {
            int k = v[1];
            pii res = query(1, 0, b.size() - 1, k);
            cout << val[lca(dp[res.first], dp[res.second])] << '\n';
        }
    }
    return;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
