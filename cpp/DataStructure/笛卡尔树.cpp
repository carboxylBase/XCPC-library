#include <bits/stdc++.h>
using namespace std;

/*
先 init 初始化
传入 vector<T> 构建小根堆笛卡尔树
*/

template<typename T>
class CartesianTree {
public:
    struct Node {
        T val;
        int ch[2];
        Node() {
            val = ch[0] = ch[1] = 0;
        }
    };
    vector<Node> nodes;

    void init(int N_) {
        nodes.clear();
        nodes.resize(N_ + 1);
    }

    void buildMinHeap(vector<T>& a) {
        init((int)a.size() - 1);
        vector<int> stk;
        stk.push_back(0);
        // a 为 1 下标数组；0 仅作栈哨兵，不比较 a[0]。
        for (int i = 1;i<a.size();i++) {
            int lst = -1;
            while (stk.size() > 1 && a[stk.back()] > a[i]) {
                lst = stk.back();
                stk.pop_back();
            }

            if (lst != -1) {
                nodes[i].ch[0] = lst;
            }
            nodes[i].val = a[i];
            nodes[stk.back()].ch[1] = i;
            stk.push_back(i);
        }
        return;
    }
};