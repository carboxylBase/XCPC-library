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
const ll MOD = 998244353;

ll qpow(ll base,ll k,ll mod) {
    if (base == 0) return 0;
    ll res = 1;
    base %= mod; base = (base + mod) % mod;
    k %= (mod - 1); k = (k + mod - 1) % (mod - 1);
    while (k) {
        if (k & 1) {
            res *= base; res %= mod;
        }
        k >>= 1;
        base *= base; base %= mod;
    }
    return res;
}

#define MMOD(x) ((x % MOD + MOD) % MOD)
struct Complex3 {
    ll x, y;
    Complex3(ll X_ = 0, ll Y_ = 0) {
        x = X_, y = Y_;
    }
    Complex3 operator+(const Complex3& z) const {
        return {MMOD(x + z.x), MMOD(y + z.y)};
    }
    Complex3 operator-(const Complex3& z) const {
        return {MMOD(x - z.x), MMOD(y - z.y)};
    }
    Complex3 operator*(const Complex3& z) const {
        ll a = x, b = y, c = z.x, d = z.y;
        return {MMOD(a*c%MOD - b*d%MOD), MMOD(a*d%MOD + b*c%MOD - b*d%MOD)};
    }
    Complex3 inv() const {
        ll a = x, b = y;
        ll n = a*a - a*b + b*b;
        n = MMOD(n);
        n = qpow(n, MOD - 2, MOD);
        return {MMOD((a - b) * n), MMOD((-b) * n)};
    }
    Complex3 operator/(const Complex3& z) const {
        return (*this) * z.inv();
    }
    bool operator==(const Complex3& z) const {
        return x == z.x && y == z.y;
    }
    Complex3 operator*(const ll& z) const {
        return Complex3(MMOD(x * z), MMOD(y * z));
    }
};
using Z = Complex3;

void solve() {
    int n, m; cin >> n >> m;
    vector<vector<Z>> mat(n + 1, vector<Z>(n + 1));
    vector<array<ll, 4>> e(m);
    for (auto&v : e) cin >> v[0] >> v[1] >> v[2] >> v[3];
    for (auto [u, v, w, t] : e) {
        Z val;
        if (t == 0) {
            val.x = w;
        } else if (t == 1) {
            val.x = w;
        } else {
            val.x = w;
        }
        mat[u][u] = mat[u][u] + val;
        mat[v][v] = mat[v][v] + val;
        mat[u][v] = mat[u][v] - val;
        mat[v][u] = mat[v][u] - val;
    }
    Z ans;
    auto cal = [&]() -> Z {
        int cnt = 0;
        for (int i = 1;i<=n-1;i++) {
            if (mat[i][i] == Z()) {
                for (int j = i + 1;j<=n-1;j++) {
                    if (mat[j][i] == Z()) {

                    } else {
                        cnt ^= 1;
                        swap(mat[i], mat[j]);
                        break;
                    }
                }
            }
            if  (mat[i][i] == Z()) {
                return Z();
            }
            for (int j = i + 1;j<=n-1;j++) {
                Z P = mat[j][i] / mat[i][i];
                for (int k = i;k<=n-1;k++) {
                    mat[j][k] = mat[j][k] - P * mat[i][k];
                }
            }
        }
        Z ans;
        ans.x = 1;
        for (int i = 1;i<=n-1;i++) {
            ans = ans * mat[i][i];
        }
        if (cnt) {
            ans = Z() - ans;
        }
        return ans;
    };
    ans = ans + cal();
    for (int i = 1;i<=n;i++) {
        for (int j = 1;j<=n;j++) {
            mat[i][j] = Z();
        }
    }
    for (auto [u, v, w, t] : e) {
        Z val;
        if (t == 0) {
            val.y = w;
        } else if (t == 1) {
            val.y = 1;
            val = val.inv();
            val = val * w;
        } else {
            val.x = w;
        }
        mat[u][u] = mat[u][u] + val;
        mat[v][v] = mat[v][v] + val;
        mat[u][v] = mat[u][v] - val;
        mat[v][u] = mat[v][u] - val;
    }
    ans = ans + cal();
    for (int i = 1;i<=n;i++) {
        for (int j = 1;j<=n;j++) {
            mat[i][j] = Z();
        }
    }
    for (auto [u, v, w, t] : e) {
        Z val;
        if (t == 0) {
            val.y = 1;
            val = val * val;
            val = val * w;
        } else if (t == 1) {
            val.y = 1;
            val = val * val;
            val = val.inv();
            val = val * w;
        } else {
            val.x = w;
        }
        mat[u][u] = mat[u][u] + val;
        mat[v][v] = mat[v][v] + val;
        mat[u][v] = mat[u][v] - val;
        mat[v][u] = mat[v][u] - val;
    }
    ans = ans + cal();
    ans.x = ans.x * qpow(3, MOD - 2, MOD) % MOD;
    cout << ans.x << endl;
    return;
}

signed main() {
#if DEBUG
    freopen("input.txt", "r", stdin);
    auto start_time = chrono::steady_clock::now();
#else
    ios::sync_with_stdio(false);
#endif
    cin.tie(nullptr);

    int t = 1;
    // cin >> t;

    while (t--) {
        solve();
    }

#if DEBUG
    auto end_time = chrono::steady_clock::now();
    auto diff = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cerr << "Time: " << diff.count() << " ms" << endl;

#endif

    return 0;
}
