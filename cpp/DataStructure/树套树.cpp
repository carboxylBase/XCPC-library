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
const ll M = 1e8 + 10;

/*
这是一个树状数组套线段树的代码
树套树的细节还是非常多的
得写的很细很细才可以
*/

int n;
namespace TWT {
    #define ls(rt) nodes[rt].ch[0]
    #define rs(rt) nodes[rt].ch[1]
    struct Node {
        int ch[2], sum;
        Node() {
            ch[0] = ch[1] = sum = 0;
        }
    } nodes[N * 4];

    int tp, rt[N];
    void init() {
        tp = 0;
        for (int i = 1;i<n+1;i++) {
            rt[i] = 0;
        }
    }

    int New() {
        nodes[++ tp] = Node();
        return tp;
    }

    void push_up(int rt) {
        nodes[rt].sum = 0;
        if (ls(rt)) {
            nodes[rt].sum += nodes[ls(rt)].sum; 
        }
        if (rs(rt)) {
            nodes[rt].sum += nodes[rs(rt)].sum;
        }
        return;
    }

    void seg_update(int& rt, int l, int r, int q, int v) {
        if (!rt) rt = New();
        if (l == r) {
            nodes[rt].sum += v;
            return;
        }
        int mid = l + r >> 1;
        if (q <= mid) seg_update(ls(rt), l, mid, q, v);
        else seg_update(rs(rt), mid + 1, r, q, v);
        push_up(rt);
        return;
    }

    int seg_query(int rt, int l, int r, int ql, int qr) {
        if (!rt) return 0;
        if (qr < l || r < ql || ql > qr) return 0;
        if (ql <= l && r <= qr) return nodes[rt].sum;
        int mid = l + r >> 1;
        return seg_query(ls(rt), l, mid, ql, qr) + seg_query(rs(rt), mid + 1, r, ql, qr);
    }

    int lowbit(int x) { return x & -x; }

    void update(int x, int v, int c) {
        while (x <= n) {
            seg_update(rt[x], 0, M, v, c);
            x += lowbit(x);
        }
        return;
    }

    int query(int x, int k) {
        if (x == 0) return 0;
        int ans = 0;
        while (x) {
            ans += seg_query(rt[x], 0, M, 0, k - 1);
            x -= lowbit(x);
        }
        return ans;
    }

    int query(int l, int r, int k) {
        return query(r, k) - query(l - 1, k);
    }

    int L[N], R[N], tpL, tpR;
    int findK(int l, int r, int k) {
        tpL = tpR = 0;
        while (r) {
            R[tpR ++] = rt[r];
            r -= lowbit(r);
        }
        l --;
        while (l) {
            L[tpL ++] = rt[l];
            l -= lowbit(l);
        }

        l = 0, r = M;
        while (l < r) {
            int mid = l + r >> 1;
            int tot = 0;
            for (int j = 0;j<tpL;j++) {
                tot -= nodes[ls(L[j])].sum;
            }
            for (int j = 0;j<tpR;j++) {
                tot += nodes[ls(R[j])].sum;
            }
            if (tot >= k) {
                for (int j = 0;j<tpL;j++) {
                    L[j] = ls(L[j]);
                }
                for (int j = 0;j<tpR;j++) {
                    R[j] = ls(R[j]);
                }
                r = mid;
            } else {
                k -= tot;
                for (int j = 0;j<tpL;j++) {
                    L[j] = rs(L[j]);
                }
                for (int j = 0;j<tpR;j++) {
                    R[j] = rs(R[j]);
                }
                l = mid + 1;
            }
        }
        return l;
    }

    int get4(int l, int r, int k) {
        k = query(l, r, k);
        if (k == 0) return -2147483647;
        return findK(l, r, k);
    }

    int get5(int l, int r, int k) {
        k = query(l, r, k + 1);
        if (k >= r - l + 1) return 2147483647;
        int ans = findK(l, r, k + 1);
        return ans;
    }
};
using namespace TWT;
