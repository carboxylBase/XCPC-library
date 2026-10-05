#include <bits/stdc++.h>

using namespace std;
using ll = long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 998244353 - 1;


// 调用 cal(a, n, p) 输出 (a^x = n) % p
// 复杂度 sqrt(p) * log
struct ExBSGS {
    static ll modpow(ll a, ll e, ll mod){
        ll r = 1;
        a %= mod;
        while(e){
            if(e & 1) r = (i128)r * a % mod;
            a = (i128)a * a % mod;
            e >>= 1;
        }
        return r;
    }
    static ll exgcd(ll a, ll b, ll &x, ll &y){
        if(b == 0){ x = 1; y = 0; return a; }
        ll x1, y1;
        ll g = exgcd(b, a % b, x1, y1);
        x = y1;
        y = x1 - a / b * y1;
        return g;
    }
    static ll invmod(ll a, ll mod){
        ll x, y;
        ll g = exgcd(a, mod, x, y);
        if(g != 1) return -1;
        x %= mod;
        if(x < 0) x += mod;
        return x;
    }
    static ll bsgs(ll a, ll b, ll mod){
        a %= mod; b %= mod;
        if(mod == 1) return 0;
        ll m = (ll)ceil(sqrt((double)mod));
        unordered_map<ll, ll> mp;
        mp.reserve(m * 2);
        ll aj = 1;
        for(ll j = 0; j < m; ++j){
            if(mp.find(aj) == mp.end()) mp[aj] = j;
            aj = (i128)aj * a % mod;
        }
        ll factor = modpow(a, m, mod);
        ll invfactor = invmod(factor, mod);
        if(invfactor == -1) return -1;
        ll cur = b % mod;
        for(ll i = 0; i <= m; ++i){
            auto it = mp.find(cur);
            if(it != mp.end()){
                return i * m + it->second;
            }
            cur = (i128)cur * invfactor % mod;
        }
        return -1;
    }
    static ll cal(ll a, ll n, ll p){
        if(p == 1) return 0;
        a %= p; n %= p;
        if(n == 1) return 0;
        ll cnt = 0;
        ll t = 1;
        ll g;
        while((g = std::gcd(a, p)) > 1){
            if(n == t) return cnt;
            if(n % g != 0) return -1;
            p /= g;
            n /= g;
            t = (i128)t * (a / g) % p;
            ++cnt;
        }
        ll invt = invmod(t, p);
        if(invt == -1) return -1;
        ll rhs = (i128)n * invt % p;
        ll res = bsgs(a, rhs, p);
        if(res == -1) return -1;
        return res + cnt;
    }
} solver_bgsg;

const int mod = 998244352;
template <typename T>
struct Mat
{
    int n, m;
    T **a;
    Mat(int _n = 0, int _m = 0) : n(_n), m(_m)
    {
        a = new T *[n];
        for (int i = 0; i < n; i++)
            a[i] = new T[m], memset(a[i], 0, sizeof(T) * m);
    }
    Mat(const Mat &B)
    {
        n = B.n, m = B.m;
        a = new T *[n];
        for (int i = 0; i < n; i++)
            a[i] = new T[m], memcpy(a[i], B.a[i], sizeof(T) * m);
    }
    ~Mat() {
        for (int i = 0; i < n; i++) delete[] a[i];
        delete[] a;
    }
    Mat &operator=(const Mat &B)
    {
        if (this == &B) return *this;
        for (int i = 0; i < n; i++) delete[] a[i];
        delete[] a;
        n = B.n, m = B.m;
        a = new T *[n];
        for (int i = 0; i < n; i++)
            a[i] = new T[m], memcpy(a[i], B.a[i], sizeof(T) * m);
        return *this;
    }
    Mat operator+(const Mat &B) const
    {
        assert(n == B.n && m == B.m);
        Mat ret(n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                ret.a[i][j] = (a[i][j] + B.a[i][j]) % mod;
        return ret;
    }
    Mat &operator+=(const Mat &B) { return *this = *this + B; }
    Mat operator*(const Mat &B) const
    {
        Mat ret(n, B.m);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < B.m; ret.a[i][j++] %= mod)
                for (int k = 0; k < m; ++k)
                    ret.a[i][j] += a[i][k] * B.a[k][j] % mod;
        return ret;
    }
    Mat &operator*=(const Mat &B) { return *this = *this * B; }
};
Mat<ll> qpow(Mat<ll> A, ll b)
{
    Mat<ll> ret(A);
    for (--b; b; b >>= 1, A *= A)
        if (b & 1)
            ret *= A;
    return ret;
}

ll exgcd(ll a, ll b, ll &x, ll &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    ll u, v;
    ll g = exgcd(b, a % b, u, v);
    x = v;
    y = u - (i128)(a / b) * v;
    return g;
}

// 求 k * x ≡ m (mod n)
// 返回最小非负解，无解返回 -1
ll solve_congruence(ll k, ll m, ll n) {
    assert(n > 0);

    k %= n;
    m %= n;
    if (k < 0) k += n;
    if (m < 0) m += n;

    ll x, y;
    ll g = exgcd(k, n, x, y);

    if (m % g != 0) return -1;

    ll period = n / g;
    ll ans = (i128)x * (m / g) % period;
    if (ans < 0) ans += period;
    return ans;
}

void solve() {
    int k; cin >> k;
    vector<ll> b(k + 1, 0);
    for (int i = 1;i<k+1;i++) {
        cin >> b[i];
    }

    ll n, m; cin >> n >> m;
    m = solver_bgsg.cal(3, m, MOD + 1);

    Mat<ll> A(k, k);
    for (int i = 0;i<k-1;i++) {
        A.a[i][i + 1] = 1;
    }
    for (int i = 0;i<k;i++) {
        A.a[k - 1][i] = b[k - i];
    }
    A = qpow(A, n - k);

    ll rhs = 0;
    rhs = A.a[k - 1][k - 1];

    ll x = solve_congruence(rhs, m, MOD);

    if (x == -1) {
        cout << x << endl;
    } else {
        x = solver_bgsg.modpow(3, x, MOD + 1);
        cout << x << endl;
    }

    return;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
}
