/**
 * G. Idiot First Search and Queries
 * * 修正点：
 * 1. 输入读取：题目中的 n 行分别对应节点 1 到 n 的子节点信息。
 * 2. 逻辑统一：一旦 Bob 离开当前子树 v 进入父节点 p，
 * p 处于"未访问"状态（因为之前是第一次从 p 下到 v，现在回溯实际上触发了 p 的重置逻辑）。
 * 这意味着进入 p 后，等同于在 p 进行一次完整的 Idiot Search。
 * 因此，我们可以把初始节点 v 的遍历也视为这个过程的第一步。
 */

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int MAXN = 300005;
const int LOGN = 20;

int n, q;
int l_child[MAXN], r_child[MAXN];
int parent[MAXN];
long long sz[MAXN];       
int tin[MAXN];            
vector<int> tour_nodes;   
int up[MAXN][LOGN];       
long long lift_cost[MAXN][LOGN]; // 表示在祖先节点 u 进行一次完整遍历所需的步数

// 构建欧拉序列 (Euler Tour) 和子树大小
void dfs_build(int u) {
    sz[u] = 1;
    tin[u] = tour_nodes.size();
    tour_nodes.push_back(u);

    if (l_child[u]) {
        dfs_build(l_child[u]);
        sz[u] += sz[l_child[u]];
        tour_nodes.push_back(u); 
    }
    
    if (r_child[u]) {
        dfs_build(r_child[u]);
        sz[u] += sz[r_child[u]];
        tour_nodes.push_back(u); 
    }
}

void solve() {
    if (!(cin >> n >> q)) return;

    // 初始化
    // 节点编号 0 到 n
    for (int i = 0; i <= n; i++) {
        l_child[i] = 0;
        r_child[i] = 0;
        parent[i] = 0;
        sz[i] = 0;
        // 清空倍增数组，防止多测污染
        for(int j=0; j<LOGN; j++) {
            up[i][j] = 0;
            lift_cost[i][j] = 0;
        }
    }
    tour_nodes.clear();
    // 预留空间防止 resize 带来的开销
    tour_nodes.reserve(2 * n + 5);

    // 输入读取修正：n 行对应节点 1 到 n
    for (int i = 1; i <= n; i++) {
        int l, r;
        cin >> l >> r;
        l_child[i] = l;
        r_child[i] = r;
        if (l) parent[l] = i;
        if (r) parent[r] = i;
    }
    
    // 节点 0 是 1 的父节点
    parent[1] = 0;

    // 从 1 开始构建，因为我们只需要在 1 的子树内模拟
    // 如果跳出 1 到达 0，游戏结束，但题目保证 k < Tv，所以不会跳出 1
    dfs_build(1);

    // 预处理倍增
    // lift_cost[u][0] 代表：在 u 节点进行一次完整遍历（并最终跳向父节点）所需的步数
    // 步数 = 2 * sz[u] - 1
    for (int i = 1; i <= n; i++) {
        up[i][0] = parent[i];
        lift_cost[i][0] = 2 * sz[i] - 1;
    }

    // 0 号节点作为哨兵，代价无穷大
    lift_cost[0][0] = 4e18; 

    for (int j = 1; j < LOGN; j++) {
        for (int i = 1; i <= n; i++) {
            int mid = up[i][j-1];
            if (mid != 0) {
                up[i][j] = up[mid][j]; // 注意：这里是接力跳
                // 逻辑解释：Bob 在 i 节点遍历完 -> 到了 parent(i) -> 在 parent(i) 遍历完 -> ...
                // 这里的倍增稍微不同于求 LCA。
                // 我们是从 i 跳出来，到了 parent(i)。
                // 然后需要在 parent(i) 进行一次完整遍历，跳到 parent(parent(i))。
                // 所以 lift_cost[i][j] 表示：从 i 能够连续跳跃 2^j 次所消耗的总步数
                // 第一次跳跃消耗：遍历 i 的代价。
                // 第二次跳跃消耗：遍历 parent(i) 的代价。
                // 
                // 修正倍增构建逻辑：
                // up[i][j-1] 是从 i 跳了 2^(j-1) 步到达的祖先 A。
                // lift_cost[i][j] = lift_cost[i][j-1] + lift_cost[A][j-1]。
                
                up[i][j] = up[mid][j-1];
                lift_cost[i][j] = lift_cost[i][j-1] + lift_cost[mid][j-1];
            } else {
                up[i][j] = 0;
                lift_cost[i][j] = 4e18; // 足够大的数
            }
        }
    }

    // 处理查询
    for (int j = 0; j < q; j++) {
        int v;
        long long k;
        cin >> v >> k;

        // 步骤 1: 检查 k 是否在当前 v 的这一轮遍历中
        // 在 v 的这一轮遍历中，最大下标是 2*sz[v] - 2 (对应最后一次回到 v)
        // 消耗 2*sz[v] - 1 步后，会跳到 parent[v]
        long long current_limit = 2 * sz[v] - 1; // 跳出的代价

        if (k < current_limit) {
            cout << tour_nodes[tin[v] + k] << "\n";
            continue;
        }

        // 步骤 2: 倍增向上跳
        // 我们需要找到一个祖先 u，使得我们消耗掉 k 的大部分，
        // 最终剩下的 k 落在 u 的某次完整遍历中。
        
        // 当前状态：我们要离开 v，消耗 current_limit，到达 parent[v]
        // 接下来就是一个标准的连跳问题：
        // 站在 v，我们要跳 1 步(完成v的遍历)，到达 parent[v]。
        // 站在 parent[v]，跳 1 步(完成parent[v]的遍历)，到达 parent[parent[v]]。
        // 所以我们其实是把起始点 v 也看作跳跃序列的第 0 个节点。
        
        k -= current_limit;
        int curr = parent[v]; // 现在位于 parent[v]，准备开始 parent[v] 的完整遍历

        // 尝试用倍增跳过祖先的完整遍历
        for (int bit = LOGN - 1; bit >= 0; bit--) {
            // 如果能跳，并且跳完之后的剩余 k 仍然 >= 0 (不对，是 k >= cost)
            // 实际上我们需要剩下的 k 小于当前所在节点的遍历代价。
            // 如果 k >= lift_cost[curr][bit]，说明我们可以完成这 2^bit 个祖先的遍历
            if (curr != 0 && up[curr][bit] != 0 && k >= lift_cost[curr][bit]) {
                k -= lift_cost[curr][bit];
                curr = up[curr][bit];
            }
        }
        
        // 此时 k 应该小于当前节点 curr 的遍历代价 (2*sz[curr]-1)
        // 答案就是 curr 这一轮遍历中的第 k 个位置
        cout << tour_nodes[tin[curr] + k] << "\n";
    }
}

int main() {
    freopen("input.txt", "r", stdin);
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int t;
    if (cin >> t) {
        while (t--) {
            solve();
        }
    }
    return 0;
}