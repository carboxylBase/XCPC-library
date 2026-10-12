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
黑科技
线性时间计算字符串最小表示.
*/

string Booth(string& s) {
    string t = s + s;
    int n = s.size();
    int i = 0, j = 1, k = 0, pos = 0;
    while (i < n && j < n && k < n) {
        if (t[i + k] == t[j + k]) {
            k ++;
        } else if (t[i + k] > t[j + k]) {
            i += k + 1;
            if (i <= j) i = j + 1;
            k = 0;
            pos = j;
        } else {
            j += k + 1;
            if (j <= i) j = i + 1;
            k = 0;
            pos = i;
        }
    }
    return t.substr(pos, n);
}