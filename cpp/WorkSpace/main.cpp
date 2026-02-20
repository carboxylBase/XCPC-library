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
const ll M = 3;

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
        nodes[rt].sum = max({nodes[rt].val, nodes[LS].sum, nodes[RS].sum});
        for (int i = 0; i < M; ++i) {
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
        nodes[rt].ch[0] = build(l, mid - 1, (dim + 1) % M);
        nodes[rt].ch[1] = build(mid + 1, r, (dim + 1) % M);
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

	ll query1(int rt, array<ll, M> v) {
		if (!rt) return 0;
		bool ok = 1;
		for (int i = 0;i<M;i++) {
			if (nodes[rt].mx[i] > v[i]) ok = 0;
		}
		if (ok) {
			return nodes[rt].sum;
		}
		ok = 0;
		for (int i = 0;i<M;i++) {
			if (nodes[rt].mn[i] > v[i]) ok = 1;
		}
		if (ok) {
			return 0;
		}
		ok = 1;
		for (int i = 0;i<M;i++) {
			if (nodes[rt].x[i] <= v[i]) {

			} else {
				ok = 0;
			}
		}
		ll res = 0;
		if (ok) {
			res = nodes[rt].val;
		}
		return max({res, query1(LS, v), query1(RS, v)});
	}

    // 核心查询逻辑
    int query(array<ll, M> v) {
		ll res = 0;
		for (int i = 0;i<20;i++) {
			res = max(res, query1(roots[i], v));
		}
		return res;
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

void solve() {
	int n; cin >> n;
	vector<array<ll, M + 1>> a(n);
	for (int i = 0;i<n;i++) {
		for (int j = 0;j<M+1;j++) {
			cin >> a[i][j];
		}
	}
	sort(a.begin(), a.end());

	ll ans = 0;

	array<ll, M> z;
	for (int i = 0;i<n;i++) {
		z[0] = a[i][1], z[1] = a[i][2], z[2] = a[i][3];
		ll res = solver.query(z) + 1;
		ans = max(ans, res);
		solver.insert(z, res);
	}

	cout << ans << endl;
	return;
}

signed main() {
#if DEBUG
	freopen("input.txt", "r", stdin);
	auto start_time = chrono::steady_clock::now();
#else
	ios::sync_with_stdio(false);
#endif
	cin.tie(nullptr);

	int t = 1;
	// cin >> t;

	while (t--) {
		solve();
	}

#if DEBUG
	auto end_time = chrono::steady_clock::now();
	auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
	cerr << "Time: " << diff.count() << " ms" << endl;

#endif

	return 0;
}
