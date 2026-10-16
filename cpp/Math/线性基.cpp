#include <bits/stdc++.h>
using ll = long long;
const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;
using namespace std;

struct LinearBasis {
    ll basis[63];
    bool flag0; // 是否存在非空子集异或为 0；输入范围 [0, 2^61)

    LinearBasis() { init(); }
    void insert(ll x) {
        if (!x) { flag0 = true; return; }
        for (int i = 60;i>=0;i--) {
            if ((1ll<<i)&x) {
                if (basis[i]) {
                    x ^= basis[i];
                } else {
                    basis[i] = x;
                    return;
                }
            }
        }
        if (!x) flag0 = true;
        return;
    }

    void preprocess() {
        for (int i = 60;i>=0;i--) {
            if (!basis[i]) continue;
            for (int j = i-1;j>=0;j--) {
                if (basis[j] && ((1ll<<j)&basis[i])) {
                    basis[i] ^= basis[j];
                }
            }
        }
        return;
    }

    void init() {
        flag0 = 0;
        memset(basis,0,sizeof basis);
        return;
    }

    // 最小非空子集异或和；没有插入任何数时返回 -1。
    ll get_mn() {
        if (flag0) return 0;
        for (int i = 0;i<=60;i++) {
            if (basis[i]) {
                return basis[i];
            }
        }
        return -1;
    }

    ll get_mx() {
        ll res = 0;
        for (int i = 60;i>=0;i--) {
            if (basis[i] && !((1ll<<i)&res)){
                res ^= basis[i];
            }
        }
        return res;
    } 
} ;