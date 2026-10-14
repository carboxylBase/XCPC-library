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

/*
先输入 s， 然后调用 build_sa 生成 sa， 再用 build_height 生成 height 数组
*/

struct SA {
    string s;
    int n;
    vector<int> sa, rk, height;

    void build_sa() {
        n = s.size();
        sa.assign(n + 1, 0);
        rk.assign(2 * n + 2, 0);
        height.assign(n + 1, 0);
        vector<int> oldrk(2 * n + 2, 0);

        for (int i = 1; i <= n; i++) {
            sa[i] = i;
            rk[i] = (unsigned char)s[i - 1] + 1;
        }

        for (int w = 1; w < n; w <<= 1) {
            sort(sa.begin() + 1, sa.end(), [&](int x, int y) {
                if (rk[x] != rk[y]) return rk[x] < rk[y];
                return rk[x + w] < rk[y + w];
            });

            oldrk = rk;
            int p = 0;
            for (int i = 1; i <= n; i++) {
                int x = sa[i], y = sa[i - 1];
                if (i == 1 || oldrk[x] != oldrk[y] ||
                    oldrk[x + w] != oldrk[y + w]) {
                    ++p;
                }
                rk[x] = p;
            }
            if (p == n) break;
        }

        // 保证 rk 为排名，也处理 n = 1 的情况
        for (int i = 1; i <= n; i++) rk[sa[i]] = i;
    }

    void build_height() {
        height.assign(n + 1, 0);
        for (int i = 1, h = 0; i <= n; i++) {
            if (rk[i] == 1) {
                h = 0;
                continue;
            }
            if (h) --h;
            int j = sa[rk[i] - 1];
            while (i + h <= n && j + h <= n &&
                   s[i + h - 1] == s[j + h - 1]) {
                ++h;
            }
            height[rk[i]] = h;
        }
    }
} solver;

int fa[N], L[N], R[N];

int F(int x) {
    if (x == fa[x]) return x;
    return fa[x] = F(fa[x]);
}

void merge(int x, int y) {
    x = F(x), y = F(y);
    L[y] += L[x], R[y] += R[x];
    fa[x] = y;
    return;
}

struct Manacher {
    ll cal(string& s) {
        string t;
        char z = '#';
        t.push_back(z); // 注意,一定填写一个 没出现过的字符
        for (auto v : s) {
            t.push_back(v); t.push_back(z);
        }
        int mid = 0, mr = 0;
        vector<int> p(t.size(), 0);
        for (int i = 1;i<t.size();i++) {
            if (mid <= i && i <= mr) {
                int j = 2 * mid - i;
                if (p[j] + i < mr) {
                    p[i] = p[j];
                    continue;
                }
                p[i] = mr - i;
                for (int j = p[i] + 1;;j++) {
                    if (i+j<t.size() && i-j>=0 && t[i+j] == t[i-j]) {

                    } else {
                        p[i] = j - 1;
                        mr = i + j - 1;
                        mid = i;
                        break;
                    }
                }
            } else {
                for (int j = 1;;j++) {
                    if (i+j<t.size() && i-j>=0 && t[i+j] == t[i-j]) {

                    } else {
                        p[i] = j - 1;
                        mr = i + j - 1;
                        mid = i;
                        break;
                    }
                }
            }
        }
        ll ans = 0;
        for (int i = 0;i<t.size();i++) {
            ll res;
            if (i & 1) {
                res = p[i] / 2 * 2 + 1;
            } else {
                res = (p[i] + 1) / 2 * 2;
            }
            ans += (res + 1) / 2;
        }
        return ans;
    }
} ma;

void solve() {
    string s; cin >> s;
    int n = s.size();

    solver.s = s;
    solver.s.push_back('#');
    reverse(s.begin(), s.end());
    solver.s += s;
    solver.build_sa();
    solver.build_height();

    struct Edge {
        int u, v, w;
        bool operator<(const Edge& A) const {
            return w < A.w;
        }
    };
    priority_queue<Edge> e;

    for (int i = 1;i+1<solver.sa.size();i++) {
        Edge z;
        z.u = solver.sa[i], z.v = solver.sa[i + 1];
        z.w = solver.height[i + 1];
        e.push(z);
    }

    for (int i = 1;i<solver.sa.size();i++) {
        fa[i] = i;
        L[i] = R[i] = 0;
        if (i <= n) L[i] = 1;
        else if (i > n + 1) R[i] = 1;
    }

    ll C = 0;

    while (!e.empty()) {
        Edge z = e.top();
        e.pop();

        int u = F(z.u), v = F(z.v), w = z.w;
        if (u == v) continue;

        C += (1ll * L[u] * R[v] + 1ll * R[u] * L[v]) * w;
        merge(u, v);
    }

    C += ma.cal(s);
    cout << C / 2 << endl;
    return;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
