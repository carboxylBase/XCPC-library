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

// max_prime_factor, 输入 n, 返回 n 的最大质因子
// 如果 n 是质数, 返回 -1

namespace primeCheck {
    std::mt19937_64 rng((unsigned)std::chrono::steady_clock::now().time_since_epoch().count());
    ll mulmod(ll a, ll b, ll mod){ return (unsigned __int128)a*b%mod; }
    ll powmod(ll a, ll d, ll mod){
        ll r=1;
        while(d){
            if(d&1) r=mulmod(r,a,mod);
            a=mulmod(a,a,mod);
            d>>=1;
        }
        return r;
    }
    bool isPrime(ll n){
        if(n<2) return false;
        for(ll p:{2,3,5,7,11,13,17,19,23,29,31,37}){ if(n%p==0) return n==p; }
        ll d=n-1, s=0;
        while((d&1)==0){ d>>=1; ++s; }
        ll bases[] = {2,325,9375,28178,450775,9780504,1795265022};
        for(ll a:bases){
            if(a%n==0) continue;
            ll x=powmod(a,d,n);
            if(x==1||x==n-1) continue;
            bool comp=true;
            for(ll r=1;r<s;++r){
                x=mulmod(x,x,n);
                if(x==n-1){ comp=false; break; }
            }
            if(comp) return false;
        }
        return true;
    }
    ll pollards_rho(ll n){
        if(n%2==0) return 2;
        if(n%3==0) return 3;
        std::uniform_int_distribution<ll> dist(2, n-2);
        while(true){
            ll c = dist(rng);
            auto f = [&](ll x){ return (mulmod(x,x,n) + c) % n; };
            ll x = dist(rng), y = x, d = 1;
            while(d==1){
                x = f(x);
                y = f(f(y));
                d = std::gcd(std::llabs(x-y), n);
            }
            if(d!=n) return d;
        }
    }
    void factor_rec(ll n, std::vector<ll>& out){
        if(n==1) return;
        if(isPrime(n)){ out.push_back(n); return; }
        ll d = pollards_rho(n);
        factor_rec(d, out);
        factor_rec(n/d, out);
    }
    ll max_prime_factor(ll n){
        if(n<=1) return -1;
        if(isPrime(n)) return -1;
        std::vector<ll> f;
        factor_rec(n,f);
        return *std::max_element(f.begin(), f.end());
    }
}
