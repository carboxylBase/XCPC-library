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

// kd-tree, 好人类智慧。。。
// kd-tree, 好好写。。。
// 我就不信 xcpc 有人敢出 kd-tree

/*
二进制分组优化的 2D-Tree
二进制优化分组是为了保持在线情形的复杂度
不需要在线维护就直接一颗树就行了
注意， nodes[0] 是哨兵节点不要用！！！
插入均摊 log^2
查询均摊 n^0.5
*/

struct KD_Tree {
    struct Node {
        int x[2];               // 坐标
        int minC[2], maxC[2];   // 子树边界矩形
        int val, sum;           // 点权值，子树总和
        int ch[2];              // 左右儿子

        void init(int _x, int _y, int _v) {
            x[0] = minC[0] = maxC[0] = _x;
            x[1] = minC[1] = maxC[1] = _y;
            val = sum = _v;
            ch[0] = ch[1] = 0;
        }
    } t[N];

    int nodeCnt = 0;
    int stk[N], tp;

    // 更新节点的统计信息
    void pushup(int p) {
        t[p].sum = t[p].val + t[t[p].ch[0]].sum + t[t[p].ch[1]].sum;
        for (int i = 0; i < 2; ++i) {
            t[p].minC[i] = t[p].maxC[i] = t[p].x[i];
            if (t[p].ch[0]) {
                t[p].minC[i] = min(t[p].minC[i], t[t[p].ch[0]].minC[i]);
                t[p].maxC[i] = max(t[p].maxC[i], t[t[p].ch[0]].maxC[i]);
            }
            if (t[p].ch[1]) {
                t[p].minC[i] = min(t[p].minC[i], t[t[p].ch[1]].minC[i]);
                t[p].maxC[i] = max(t[p].maxC[i], t[t[p].ch[1]].maxC[i]);
            }
        }
    }

    // 递归构建平衡 KD-Tree
    int build(int l, int r, int dim) {
        if (l > r) return 0;
        int mid = (l + r) >> 1;
        nth_element(stk + l, stk + mid, stk + r + 1, [&](int a, int b) {
            return t[a].x[dim] < t[b].x[dim];
        });
        int p = stk[mid];
        t[p].ch[0] = build(l, mid - 1, dim ^ 1);
        t[p].ch[1] = build(mid + 1, r, dim ^ 1);
        pushup(p);
        return p;
    }

    // 展平子树
    void flatten(int p) {
        if (!p) return;
        stk[++tp] = p;
        flatten(t[p].ch[0]);
        flatten(t[p].ch[1]);
    }

    // 核心查询逻辑
    int query(int p, int x1, int y1, int x2, int y2) {
        if (!p) return 0;
        // 1. 完全包含
        if (t[p].minC[0] >= x1 && t[p].maxC[0] <= x2 && t[p].minC[1] >= y1 && t[p].maxC[1] <= y2)
            return t[p].sum;
        // 2. 完全不相交
        if (t[p].minC[0] > x2 || t[p].maxC[0] < x1 || t[p].minC[1] > y2 || t[p].maxC[1] < y1)
            return 0;
        // 3. 部分重叠
        int res = 0;
        if (t[p].x[0] >= x1 && t[p].x[0] <= x2 && t[p].x[1] >= y1 && t[p].x[1] <= y2)
            res += t[p].val;
        return res + query(t[p].ch[0], x1, y1, x2, y2) + query(t[p].ch[1], x1, y1, x2, y2);
    }

    // 二进制分组管理器
    int roots[25], treeSize[25];
    void insert(int x, int y, int v) {
        int cur = ++nodeCnt;
        t[cur].init(x, y, v);
        
        tp = 0;
        stk[++tp] = cur;
        
        // 类似二进制加法：1 + 1 = 10 (进位合并)
        for (int i = 0; i < 20; ++i) {
            if (!roots[i]) {
                roots[i] = build(1, tp, 0);
                treeSize[i] = tp;
                return;
            }
            flatten(roots[i]);
            roots[i] = 0; // 清空当前层，准备合并到下一层
            treeSize[i] = 0;
        }
    }
} solver;

// 这种写法， 就很优美!!
const int M = 2;
struct KD_Tree {
    struct Node {
        ll x[M], mn[M], mx[M], ch[M], val, sum;

        void init(array<ll, M> y, ll _v) {
			for (int i = 0;i<M;i++) {
				x[i] = mn[i] = mx[i] = y[i];
			}
            val = sum = _v;
            ch[0] = ch[1] = 0;
        }
    } nodes[N];

    int nodeCnt = 0;
    int stk[N], tp;

	#define LS (nodes[rt].ch[0])
	#define RS (nodes[rt].ch[1])

    // 更新节点的统计信息
    void pushup(int rt) {
        nodes[rt].sum = nodes[rt].val + nodes[LS].sum + nodes[RS].sum;
        for (int i = 0; i < 2; ++i) {
            nodes[rt].mn[i] = nodes[rt].mx[i] = nodes[rt].x[i];
            if (LS) {
                nodes[rt].mn[i] = min(nodes[rt].mn[i], nodes[LS].mn[i]);
                nodes[rt].mx[i] = max(nodes[rt].mx[i], nodes[LS].mx[i]);
            }
            if (RS) {
                nodes[rt].mn[i] = min(nodes[rt].mn[i], nodes[RS].mn[i]);
                nodes[rt].mx[i] = max(nodes[rt].mx[i], nodes[RS].mx[i]);
            }
        }
    }

    // 递归构建平衡 KD-Tree
    int build(int l, int r, int dim) {
        if (l > r) return 0;
        int mid = (l + r) >> 1;
        nth_element(stk + l, stk + mid, stk + r + 1, [&](int a, int b) {
            return nodes[a].x[dim] < nodes[b].x[dim];
        });
        int rt = stk[mid];
        nodes[rt].ch[0] = build(l, mid - 1, dim ^ 1);
        nodes[rt].ch[1] = build(mid + 1, r, dim ^ 1);
        pushup(rt);
        return rt;
    }

    // 展平子树
    void flatten(int rt) {
        if (!rt) return;
        stk[++tp] = rt;
        flatten(nodes[rt].ch[0]);
        flatten(nodes[rt].ch[1]);
    }

	ll get_mn(int rt, array<ll, M> v) {
		ll res = 0;
		if (nodes[rt].mx[0] < v[0]) {
			res += v[0] - nodes[rt].mx[0];
		} else if (nodes[rt].mn[0] > v[0]) {
			res += nodes[rt].mn[0] - v[0];
		}
		if (nodes[rt].mx[1] < v[1]) {
			res += v[1] - nodes[rt].mx[1];
		} else if (nodes[rt].mn[1] > v[1]) {
			res += nodes[rt].mn[1] - v[1];
		}
		return res;
	}

	ll get_mx(int rt, array<ll, M> v) {
		ll res = 0;
		res += max(abs(nodes[rt].mx[0] - v[0]), abs(nodes[rt].mn[0] - v[0]));
		res += max(abs(nodes[rt].mx[1] - v[1]), abs(nodes[rt].mn[1] - v[1]));
		return res;
	}

	void queryMn(int rt, array<ll, M> v, ll &mn) {
		if (!rt) return;
        if (nodes[rt].x[0] == v[0] && nodes[rt].x[1] == v[1]) {

		} else {
			mn = min(mn, abs(v[0] - nodes[rt].x[0]) + abs(v[1] - nodes[rt].x[1]));
		}
		ll lhs = INF, rhs = INF;
		if (LS) lhs = get_mn(LS, v);
		if (RS) rhs = get_mn(RS, v);
		if (lhs > rhs) {
			if (rhs < mn) queryMn(RS, v, mn);
			if (lhs < mn) queryMn(LS, v, mn); 
		} else {
			if (lhs < mn) queryMn(LS, v, mn); 
			if (rhs < mn) queryMn(RS, v, mn);
		}
		return;
	}

	void queryMx(int rt, array<ll, M> v, ll &mx) {
		if (!rt) return;
        mx = max(mx, abs(v[0] - nodes[rt].x[0]) + abs(v[1] - nodes[rt].x[1]));
		ll lhs = -INF, rhs = -INF;
		if (LS) lhs = get_mx(LS, v);
		if (RS) rhs = get_mx(RS, v);
		if (lhs < rhs) {
			if (rhs > mx) queryMx(RS, v, mx);
			if (lhs > mx) queryMx(LS, v, mx); 
		} else {
			if (lhs > mx) queryMx(LS, v, mx); 
			if (rhs > mx) queryMx(RS, v, mx);
		}
	}

    // 核心查询逻辑
    ll query(array<ll, M> v) {
		ll mn = INF, mx = -INF;
		queryMn(roots[0], v, mn);
		queryMx(roots[0], v, mx);
		return mx - mn;
    }

    // 二进制分组管理器
    int roots[25], treeSize[25];
    void insert(array<ll, M> y, ll v) {
        int cur = ++nodeCnt;
        nodes[cur].init(y, v);
        
        tp = 0;
        stk[++tp] = cur;
        
        // 类似二进制加法：1 + 1 = 10 (进位合并)
        for (int i = 0; i < 20; ++i) {
            if (!roots[i]) {
                roots[i] = build(1, tp, 0);
                treeSize[i] = tp;
                return;
            }
            flatten(roots[i]);
            roots[i] = 0; // 清空当前层，准备合并到下一层
            treeSize[i] = 0;
        }
    }
} solver;