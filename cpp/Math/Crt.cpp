/*
功能简短说明：
1. 输入：两个 `(r,m)` 分别表示 `x ≡ r (mod m)`；建议 `m>0` 且 `0≤r<m`。
2. 支持 **模不互质** 的情况（会处理公约数）。
3. 输出：若有解返回 `(r,M)`，表示所有解为 `x = r + k*M`（`k∈Z`），且 `0 ≤ r < M`；若无解返回 `(-1,-1)`。
4. 合并规则：如果方程组相容（即 `r2-r1` 能被 `g = gcd(m1,m2)` 整除），则合并得到模 `M = m1 * (m2 / g)`，并计算最小非负解 `r`。
5. 算法要点：用扩展欧几里得求 `g` 和系数 `x,y`，检查 `(r2-r1)%g==0`，计算位移 `t = (r2-r1)/g * x (mod m2/g)`，再算出 `r = (r1 + m1*t) mod M`。实现中用 `__int128` 防止乘法溢出。
   示例：`(2,3)` 与 `(3,5)` 合并得 `(8,15)`，因为解集是 `x ≡ 8 (mod 15)`。
*/
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using i128 = __int128;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
ll exgcd(ll a,ll b,ll &x,ll &y){
    if(b==0){x=1;y=0;return a;}
    ll x1,y1; ll g=exgcd(b,a%b,x1,y1);
    x=y1; y=x1-(a/b)*y1; return g;
}
pll crt(pll A,pll B){
    ll r1=A.first,m1=A.second,r2=B.first,m2=B.second;
    if(m1<=0||m2<=0) return {-1,-1};
    r1%=m1; if(r1<0) r1+=m1;
    r2%=m2; if(r2<0) r2+=m2;
    ll x,y; ll g=exgcd(m1,m2,x,y);
    ll d=r2-r1; if(d%g!=0) return {-1,-1};
    i128 t=(i128)(d/g)*(i128)x;
    i128 mod2=m2/g; t%=mod2; if(t<0) t+=mod2;
    i128 M=(i128)m1*mod2;
    i128 res=((i128)r1 + (i128)m1 * t) % M; if(res<0) res+=M;
    return {(ll)res,(ll)M};
}
pll crt_many(const vector<pll>& v){
    if(v.empty()) return {0,1};
    pll cur = v[0];
    for(size_t i=1;i<v.size();++i){
        cur = crt(cur,v[i]);
        if(cur.first==-1) return cur;
    }
    return cur;
}
ll inv_mod(ll a, ll b){
    ll x, y;
    ll g = exgcd(a, b, x, y);
    if(g != 1) return -1;
    x %= b;
    if(x < 0) x += b;
    return x;
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n; if(!(cin>>n)) return 0;
    vector<pll> v; v.reserve(n);
    for(int i=0;i<n;++i){ ll r,m; cin>>r>>m; v.emplace_back(r,m); }
    auto ans = crt_many(v);
    cout<<ans.first<<" "<<ans.second<<"\n";
    return 0;
}


using i64 = long long;
using i128 = __int128_t;

i64 power(i64 a, i64 b, i64 p) {
    i64 res = 1 % p;
    for (a %= p; b; b >>= 1, a = i128(a) * a % p) {
        if (b & 1) {
            res = i128(res) * a % p;
        }
    }
    return res;
}

// a * x ≡ b (mod p)，p > 0
// 返回最小非负解，无解返回 -1
i64 linear(i64 a, i64 b, i64 p) {
    a %= p;
    b %= p;
    if (a < 0) a += p;
    if (b < 0) b += p;

    i64 x, y;
    i64 g = exgcd(a, p, x, y);
    if (b % g) return -1;

    p /= g;
    x = i128(x) * (b / g) % p;
    if (x < 0) x += p;
    return x;
}

// a^x ≡ b (mod p)，要求 a,b >= 0，p > 1，gcd(a,p) = 1
// 返回最小非负解，无解返回 -1
i64 bsgs(i64 a, i64 b, i64 p) {
    a %= p;
    b %= p;
    if (b == 1) return 0;

    i64 s = std::sqrt((long double)p) + 1;
    std::unordered_map<i64, i64> mp;
    mp.reserve(2 * s);

    i64 v = 1;
    for (i64 j = 0; j < s; j++) {
        mp.emplace(v, j); // 重复时保留最小 j
        v = i128(v) * a % p;
    }

    i64 step = linear(v, 1, p); // a^s 的逆元
    for (i64 i = 0; i < s; i++) {
        auto it = mp.find(b);
        if (it != mp.end()) {
            return i * s + it->second;
        }
        b = i128(b) * step % p;
    }
    return -1;
}