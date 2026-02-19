#include <iostream>
#include <algorithm>

using namespace std;

const int MAXN = 200005; // 操作数 2e5
const int INF = 2e9;

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
} t[MAXN];

int nodeCount = 0;
int tempIndices[MAXN], ptr;

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
    nth_element(tempIndices + l, tempIndices + mid, tempIndices + r + 1, [&](int a, int b) {
        return t[a].x[dim] < t[b].x[dim];
    });
    int p = tempIndices[mid];
    t[p].ch[0] = build(l, mid - 1, dim ^ 1);
    t[p].ch[1] = build(mid + 1, r, dim ^ 1);
    pushup(p);
    return p;
}

// 展平子树
void flatten(int p) {
    if (!p) return;
    tempIndices[++ptr] = p;
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
    int cur = ++nodeCount;
    t[cur].init(x, y, v);
    
    ptr = 0;
    tempIndices[++ptr] = cur;
    
    // 类似二进制加法：1 + 1 = 10 (进位合并)
    for (int i = 0; i < 20; ++i) {
        if (!roots[i]) {
            roots[i] = build(1, ptr, 0);
            treeSize[i] = ptr;
            return;
        }
        flatten(roots[i]);
        roots[i] = 0; // 清空当前层，准备合并到下一层
        treeSize[i] = 0;
    }
}

int main() {
	freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(false); cin.tie(0);
    int n_useless; cin >> n_useless; // 棋盘大小 N 对 KD-Tree 没用
    int lastAns = 0, op;

    while (cin >> op && op != 3) {
        if (op == 1) {
            int x, y, a; cin >> x >> y >> a;
            insert(x ^ lastAns, y ^ lastAns, a ^ lastAns);
        } else if (op == 2) {
            int x1, y1, x2, y2; cin >> x1 >> y1 >> x2 >> y2;
            x1 ^= lastAns; y1 ^= lastAns; x2 ^= lastAns; y2 ^= lastAns;
            // 题目保证 x1<=x2, y1<=y2，但不放心可以加个 swap
            if (x1 > x2) swap(x1, x2);
            if (y1 > y2) swap(y1, y2);
            lastAns = 0;
            for (int i = 0; i < 20; ++i) 
                if (roots[i]) lastAns += query(roots[i], x1, y1, x2, y2);
            cout << lastAns << "\n";
        }
    }
    return 0;
}