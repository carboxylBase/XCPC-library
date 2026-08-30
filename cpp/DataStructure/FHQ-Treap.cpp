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
关于 FHQ-Treap 的易错点
注意懒标记的作用顺序
注意 unite 之后需要重新指定 rt 的值
*/

/*
以下是文艺平衡树模版, 功能是对区间进行若干次翻转, 最终输出反转后结果
这是根据 siz 进行反转的示例
大概的思想是把 1-(l-1) 的部分切开, 然后切 l - r 的, 然后懒标记
*/
struct FHQ_Treap {
    int ch[N][2], val[N], siz[N], rnd[N];
    bool tag[N];
    int tot, rt;

    FHQ_Treap() {
        tot = 0;
        rt = 0;
        memset(ch, 0, sizeof ch);
        memset(siz, 0, sizeof siz);
        memset(tag, 0, sizeof tag);
        srand((unsigned)time(NULL));
    }

    int New(int v) {
        ++tot;
        val[tot] = v;
        siz[tot] = 1;
        rnd[tot] = rand();
        tag[tot] = false;
        ch[tot][0] = ch[tot][1] = 0;
        return tot;
    }

    void push_up(int u) {
        siz[u] = siz[ch[u][0]] + siz[ch[u][1]] + 1;
    }

    void push_down(int u) {
        if (!u || !tag[u]) return;
        tag[u] = false;
        swap(ch[u][0], ch[u][1]);
        if (ch[u][0]) tag[ch[u][0]] ^= 1;
        if (ch[u][1]) tag[ch[u][1]] ^= 1;
    }

    void split(int u, int k, int &x, int &y) {
        if (!u) { x = y = 0; return; }
        push_down(u);
        if (siz[ch[u][0]] < k) {
            x = u;
            split(ch[u][1], k - siz[ch[u][0]] - 1, ch[u][1], y);
        } else {
            y = u;
            split(ch[u][0], k, x, ch[u][0]);
        }
        push_up(u);
    }

    int unite(int x, int y) {
        if (!x || !y) return x + y;
        if (rnd[x] < rnd[y]) {
            push_down(x);
            ch[x][1] = unite(ch[x][1], y);
            push_up(x);
            return x;
        } else {
            push_down(y);
            ch[y][0] = unite(x, ch[y][0]);
            push_up(y);
            return y;
        }
    }

    void insert(int v) {
        int a, b;
        split(rt, v, a, b);
        rt = unite(unite(a, New(v)), b);
    }

    void reverse(int l, int r) {
        int a, b, c;
        split(rt, l - 1, a, b);
        split(b, r - l + 1, b, c);
        tag[b] ^= 1;
        rt = unite(a, unite(b, c));
    }

    void dfs(int u) {
        if (!u) return;
        if (tag[u]) push_down(u);
        dfs(ch[u][0]);
        printf("%d ", val[u]);
        dfs(ch[u][1]);
    }

    void print() { dfs(rt); }
} solver;

/*
这个是按照值域分解的版本, 同时具有给指定值域的数加上一定值的功能
*/

struct FHQ_Treap {
    int ch[N][2], rnd[N], idx[N];
    int tot, rt;
    ll val[N], laz[N];

    FHQ_Treap() {
        tot = 0;
        rt = 0;
        memset(ch, 0, sizeof ch);
        memset(laz, 0, sizeof laz);
        srand((unsigned)time(NULL));
    }

    void init() {
        tot = rt = 0;
    }

    int New(int v, int id) {
        ++tot;
        val[tot] = v;
        rnd[tot] = rand();
        laz[tot] = 0;
        idx[tot] = id;
        ch[tot][0] = ch[tot][1] = 0;
        return tot;
    }

    void push_down(int u) {
        if (!u || !laz[u]) return;
        for (int i = 0;i<2;i++) {
            if (!ch[u][i]) continue;
            val[ch[u][i]] += laz[u];
            laz[ch[u][i]] += laz[u];
        }
        laz[u] = 0;
    }

    void split(int u, ll k, int &x, int &y) {
        if (!u) { x = y = 0; return; }
        push_down(u);
        if (val[u] <= k) {
            x = u;
            split(ch[u][1], k, ch[u][1], y);
        } else {
            y = u;
            split(ch[u][0], k, x, ch[u][0]);
        }
    }

    int unite(int x, int y) {
        if (!x || !y) return x + y;
        if (rnd[x] < rnd[y]) {
            push_down(x);
            ch[x][1] = unite(ch[x][1], y);
            return x;
        } else {
            push_down(y);
            ch[y][0] = unite(x, ch[y][0]);
            return y;
        }
    }

    void insert(int v, int id) {
        int a, b;
        split(rt, v, a, b);
        rt = unite(unite(a, New(v, id)), b);
    }

    void add(ll x) {
        int a, b;
        split(rt, x - 1, a, b);
        if (b) {
            val[b] += x;
            laz[b] += x;
        }
        rt = unite(a, b);
    }
} solver;