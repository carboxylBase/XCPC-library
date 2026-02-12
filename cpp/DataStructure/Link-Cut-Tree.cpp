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
const ll MOD = 51061;

/*
通俗 LCT 模版
*/
struct LCT {
    int n;
    int f[N], st[N];
    array<int, 2> c[N];
    char rev[N];

    void init(int n_) {
        n = n_;
        for (int i = 1;i<=n;i++) {
            f[i] = st[i] = c[i][0] = c[i][1] = 0;
            rev[i] = 0;
        }
    }

    bool nroot(int x) {
        int fa = f[x];
        return c[fa][0] == x || c[fa][1] == x;
    }

    void pushup(int x) {
        // 维护路径信息的核心操作
    }

    void pushr(int x) {
        swap(c[x][0], c[x][1]);
        rev[x] ^= 1;
    }

    void pushdown(int x) {
        if (rev[x]) {
            if (c[x][0]) {
                pushr(c[x][0]);
            }
            if (c[x][1]) {
                pushr(c[x][1]);
            }
            rev[x] = 0;
        }
    }

    void rotate(int x) {
        int y = f[x];
        int z = f[y];
        int k = (c[y][1] == x);
        int w = c[x][!k];

        if (nroot(y)) {
            if (c[z][1] == y) {
                c[z][1] = x;
            } else {
                c[z][0] = x;
            }
        }

        c[x][!k] = y;
        c[y][k] = w;

        if (w) {
            f[w] = y;
        }

        f[y] = x;
        f[x] = z;

        pushup(y);
    }

    void splay(int x) {
        int y = x;
        int top = 0;
        st[++top] = y;
        while (nroot(y)) {
            y = f[y];
            st[++top] = y;
        }
        while (top) {
            pushdown(st[top--]);
        }
        while (nroot(x)) {
            y = f[x];
            int z = f[y];
            bool zigzag = (c[y][0] == x) ^ (c[z][0] == y);
            if (nroot(y)) {
                if (zigzag) {
                    rotate(x);
                } else {
                    rotate(y);
                }
            }
            rotate(x);
        }
        pushup(x);
    }

    void access(int x) {
        int y = 0;
        while (x) {
            splay(x);
            c[x][1] = y;
            pushup(x);
            y = x;
            x = f[x];
        }
    }

    void makeroot(int x) {
        access(x);
        splay(x);
        pushr(x);
    }

    int findroot(int x) {
        access(x);
        splay(x);
        while (c[x][0]) {
            pushdown(x);
            x = c[x][0];
        }
        splay(x);
        return x;
    }

    void split(int x, int y) {
        makeroot(x);
        access(y);
        splay(y);
    }

    void link(int x, int y) {
        makeroot(x);
        if (findroot(y) != x) {
            f[x] = y;
        }
    }

    void cut(int x, int y) {
        makeroot(x);
        if (findroot(y) == x && f[y] == x && c[y][0] == 0) {
            f[y] = 0;
            c[x][1] = 0;
            pushup(x);
        }
    }

    void setval(int x, ll v) {
        splay(x);
        // sum[x] = val[x] = v;
        pushup(x);
    }

    int query(int x, int y) {
        split(x, y);
        // return sum[y];
    }
} solver;

/*
维护子树信息特供版

*/
struct LCT {
    int n;
    int f[N], st[N];
    array<int, 2> c[N];
    char rev[N];

    int s[N], si[N];

    bool nroot(int x) {
        int fa = f[x];
        return c[fa][0] == x || c[fa][1] == x;
    }

    void pushup(int x) {
        // 维护路径信息的核心操作
        // 下面是维护两种常见信息 示例
        s[x] = s[c[x][0]] + s[c[x][1]] + si[x] + 1;
        // s[x] = s[c[x][0]] + s[c[x][1]] + si[x] + v[x];
    }

    void pushr(int x) {
        swap(c[x][0], c[x][1]);
        rev[x] ^= 1;
    }

    void pushdown(int x) {
        if (rev[x]) {
            if (c[x][0]) {
                pushr(c[x][0]);
            }
            if (c[x][1]) {
                pushr(c[x][1]);
            }
            rev[x] = 0;
        }
    }

    void rotate(int x) {
        int y = f[x];
        int z = f[y];
        int k = (c[y][1] == x);
        int w = c[x][!k];

        if (nroot(y)) {
            if (c[z][1] == y) {
                c[z][1] = x;
            } else {
                c[z][0] = x;
            }
        }

        c[x][!k] = y;
        c[y][k] = w;

        if (w) {
            f[w] = y;
        }

        f[y] = x;
        f[x] = z;

        pushup(y);
    }

    void splay(int x) {
        int y = x;
        int top = 0;
        st[++top] = y;
        while (nroot(y)) {
            y = f[y];
            st[++top] = y;
        }
        while (top) {
            pushdown(st[top--]);
        }
        while (nroot(x)) {
            y = f[x];
            int z = f[y];
            bool zigzag = (c[y][0] == x) ^ (c[z][0] == y);
            if (nroot(y)) {
                if (zigzag) {
                    rotate(x);
                } else {
                    rotate(y);
                }
            }
            rotate(x);
        }
        pushup(x);
    }

    void access(int x) {
        for (int y = 0; x; y = x, x = f[x]) {
            splay(x);
            si[x] += s[c[x][1]];
            si[x] -= s[c[x][1] = y];
            pushup(x);
        }
    }

    void makeroot(int x) {
        access(x);
        splay(x);
        pushr(x);
    }

    int findroot(int x) {
        access(x);
        splay(x);
        while (c[x][0]) {
            pushdown(x);
            x = c[x][0];
        }
        splay(x);
        return x;
    }

    void split(int x, int y) {
        makeroot(x);
        access(y);
        splay(y);
    }

    void link(int x, int y) {
        makeroot(x);
        if (findroot(y) != x) {
            makeroot(y);
            f[x] = y;
            si[y] += s[x];
            pushup(y);
        }
    }

    void cut(int x, int y) {
        makeroot(x);
        if (findroot(y) == x && f[y] == x && c[y][0] == 0) {
            f[y] = 0;
            c[x][1] = 0;
            pushup(x);
        }
    }
} solver;

/*
这个 LCT 实现区间加乘
*/
struct LCT {
    int n;
    int f[N], val[N], st[N];
    array<int, 2> c[N];
    char rev[N];
    ll sum[N], sz[N], a[N], lp[N], lm[N];

    void init(int n_) {
        n = n_;
        for (int i = 1;i<=n;i++) {
            lp[i] = 0;
            lm[i] = 1;
            sum[i] = sz[i] = a[i] = 0;
            f[i] = val[i] = st[i] = c[i][0] = c[i][1] = 0;
            rev[i] = 0;
        }
    }

    bool nroot(int x) {
        int fa = f[x];
        return c[fa][0] == x || c[fa][1] == x;
    }

    void pushup(int x) {
        // 维护路径信息的核心操作
        sum[x] = (sum[c[x][0]] + sum[c[x][1]] + a[x]) % MOD;
        sz[x] = (sz[c[x][0]] + sz[c[x][1]] + 1) % MOD;
    }

    void pushr(int x) {
        swap(c[x][0], c[x][1]);
        rev[x] ^= 1;
    }

    void pushm(int x, ll c) {
        sum[x] = sum[x] * c % MOD;
        a[x] = a[x] * c % MOD;
        lm[x] = lm[x] * c % MOD;
        lp[x] = lp[x] * c % MOD;
    }

    void pusha(int x, ll c) {
        sum[x] = (sum[x] + sz[x] * c % MOD) % MOD;
        a[x] = (a[x] + c) % MOD;
        lp[x] = (lp[x] + c) % MOD;
    }

    void pushdown(int x) {
        if (lm[x] != 1) {
            pushm(c[x][0], lm[x]);
            pushm(c[x][1], lm[x]);
            lm[x] = 1;
        }
        if (lp[x]) {
            pusha(c[x][0], lp[x]);
            pusha(c[x][1], lp[x]);
            lp[x] = 0;
        }
        if (rev[x]) {
            if (c[x][0]) {
                pushr(c[x][0]);
            }
            if (c[x][1]) {
                pushr(c[x][1]);
            }
            rev[x] = 0;
        }
    }

    void rotate(int x) {
        int y = f[x];
        int z = f[y];
        int k = (c[y][1] == x);
        int w = c[x][!k];

        if (nroot(y)) {
            if (c[z][1] == y) {
                c[z][1] = x;
            } else {
                c[z][0] = x;
            }
        }

        c[x][!k] = y;
        c[y][k] = w;

        if (w) {
            f[w] = y;
        }

        f[y] = x;
        f[x] = z;

        pushup(y);
    }

    void splay(int x) {
        int y = x;
        int top = 0;
        st[++top] = y;
        while (nroot(y)) {
            y = f[y];
            st[++top] = y;
        }
        while (top) {
            pushdown(st[top--]);
        }
        while (nroot(x)) {
            y = f[x];
            int z = f[y];
            bool zigzag = (c[y][0] == x) ^ (c[z][0] == y);
            if (nroot(y)) {
                if (zigzag) {
                    rotate(x);
                } else {
                    rotate(y);
                }
            }
            rotate(x);
        }
        pushup(x);
    }

    void access(int x) {
        int y = 0;
        while (x) {
            splay(x);
            c[x][1] = y;
            pushup(x);
            y = x;
            x = f[x];
        }
    }

    void makeroot(int x) {
        access(x);
        splay(x);
        pushr(x);
    }

    int findroot(int x) {
        access(x);
        splay(x);
        while (c[x][0]) {
            pushdown(x);
            x = c[x][0];
        }
        splay(x);
        return x;
    }

    void split(int x, int y) {
        makeroot(x);
        access(y);
        splay(y);
    }

    void link(int x, int y) {
        makeroot(x);
        if (findroot(y) != x) {
            f[x] = y;
        }
    }

    void cut(int x, int y) {
        makeroot(x);
        if (findroot(y) == x && f[y] == x && c[y][0] == 0) {
            f[y] = 0;
            c[x][1] = 0;
            pushup(x);
        }
    }

    void pls(int x, int y, ll v) {
        split(x, y);
        pusha(y, v);
    }

    void mul(int x, int y, ll v) {
        split(x, y);
        pushm(y, v);
    }

    void setval(int x, ll v) {
        splay(x);
        sum[x] = val[x] = v;
        pushup(x);
    }

    int query(int x, int y) {
        split(x, y);
        return sum[y];
    }
} solver;