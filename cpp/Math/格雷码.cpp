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

vector<int> ans;
// 生成长度为 L 的完整格雷码，从 head 开始
// 返回局部终点；L 必须是 2 的幂
int emitGray(int L, int head, int offset) {
    for (int i = 0;i<L;i++) {
        int gray = i ^ (i >> 1);
        ans.push_back(offset + (gray ^ head));
    }
    return head ^ (L / 2);
}