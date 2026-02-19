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
二进制分组优化的 2D-Tree
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

void solve() {
    int n; cin >> n;
    int lst = 0, opt;
    while (cin >> opt && opt != 3) {
        if (opt == 1) {
            int x, y, v; cin >> x >> y >> v;
            solver.insert(x ^ lst, y ^ lst, v ^ lst);
        } else {
            int x1, y1, x2, y2; cin >> x1 >> y1 >> x2 >> y2;
            x1 ^= lst, y1 ^= lst, x2 ^= lst, y2 ^= lst;
            lst = 0;
            for (int i = 0;i<20;i++) {
                if (solver.roots[i]) {
                    lst += solver.query(solver.roots[i], x1, y1, x2, y2);
                }
            }
            cout << lst << endl;
        }
    }
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
