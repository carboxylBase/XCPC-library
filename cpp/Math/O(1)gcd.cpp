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

// 先调用 init 函数！！！
// 复杂度是线性的！！！

namespace quickGcd {
    const int M = 1000000;
    const int T = 1000;
    
    int pre[T + 1][T + 1];
    int fac[M + 1][3];
    int pri[M / 10], tot;
    bool not_p[M + 1];

    // 预处理函数：在线性筛的基础上进行因子分解和GCD表制作
    void init() {
        for (int i = 0; i <= T; ++i) {
            pre[i][0] = pre[0][i] = i;
            for (int j = 1; j <= i; ++j)
                pre[i][j] = pre[j][i] = pre[j][i % j];
        }

        fac[1][0] = fac[1][1] = fac[1][2] = 1;
        for (int i = 2; i <= M; ++i) {
            if (!not_p[i]) {
                pri[++tot] = i;
                fac[i][0] = fac[i][1] = 1;
                fac[i][2] = i;
            }
            for (int j = 1; j <= tot && i * pri[j] <= M; ++j) {
                int tmp = i * pri[j];
                not_p[tmp] = true;
                fac[tmp][0] = fac[i][0] * pri[j];
                fac[tmp][1] = fac[i][1];
                fac[tmp][2] = fac[i][2];
                
                // 保持 fac[0] <= fac[1] <= fac[2]
                if (fac[tmp][0] > fac[tmp][1]) std::swap(fac[tmp][0], fac[tmp][1]);
                if (fac[tmp][1] > fac[tmp][2]) std::swap(fac[tmp][1], fac[tmp][2]);
                if (fac[tmp][0] > fac[tmp][1]) std::swap(fac[tmp][0], fac[tmp][1]);

                if (i % pri[j] == 0) break;
            }
        }
    }

    // O(1) 查询接口
    inline int gcd(int a, int b) {
        int res = 1;
        for (int i = 0; i < 3; ++i) {
            int f = fac[a][i];
            int cur;
            if (f <= T) {
                cur = pre[f][b % f];
            } else {
                cur = (b % f == 0) ? f : 1;
            }
            b /= cur;
            res *= cur;
        }
        return res;
    }
};
