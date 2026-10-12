#include "all.hpp"
using namespace std;
#define DEBUG 1
using ll = long long;
using pii = pair<ll,ll>;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;

/*
Warning: 
1. 统一使用 -1 表示无解
2. 假定均为小写字母

Usage:
1. init()
2. insert() 所有字符串
3. process() 预处理
*/

namespace ACAM {
    struct Node {
        int ch[26], idx, fail, out;
        Node() {
            for (int i = 0;i<26;i++) ch[i] = -1;
            idx = fail = out = -1;
        }
    } nodes[N];

    int rt, cnt_node;

    int newNode() {
        nodes[cnt_node] = Node();
        return cnt_node ++;
    }

    void init() {
        cnt_node = 0;
        rt = newNode();
        return;
    }

    // 这里默认使用 ACAM::rt, dep = 0
    void insert(int rt, int dep, string& s, int idx) {
        if (s.size() == dep) {
            nodes[rt].idx = idx;
            return;
        }
        if (nodes[rt].ch[s[dep] - 'a'] == -1) nodes[rt].ch[s[dep] - 'a'] = newNode();
        insert(nodes[rt].ch[s[dep] - 'a'], dep + 1, s, idx);
        return;
    }

    void process() {
        queue<int> q;
        
        for (int i = 0;i<26;i++) {
            int v = nodes[rt].ch[i];
            if (v != -1) {
                q.push(v);
                nodes[v].fail = rt;
            } else {
                nodes[rt].ch[i] = rt;
            }
        }

        while (!q.empty()) {
            int u = q.front(); 
            q.pop();

            for (int i = 0;i<26;i++) {
                int v = nodes[u].ch[i];
                if (v == -1) {
                    nodes[u].ch[i] = nodes[nodes[u].fail].ch[i];
                } else {
                    q.push(v);
                    nodes[v].fail = nodes[nodes[u].fail].ch[i];
                    if (nodes[nodes[v].fail].idx != -1) {
                        nodes[v].out = nodes[v].fail;
                    } else {
                        nodes[v].out = nodes[nodes[v].fail].out;
                    }
                }
            }
        }

        return;
    }
};

