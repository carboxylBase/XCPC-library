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
const int M = 17;

// 值得注意的是， 这个输出的结果似乎没有取模
template <typename T>
struct FWTSolver {
    enum Type { OR = 0, AND = 1, XOR = 2 };

    // 核心变换
    static void transform(vector<T>& a, Type type, bool inv) {
        int n = a.size();
        for (int h = 1; h < n; h <<= 1) {
            for (int j = 0; j < n; j += (h << 1)) {
                for (int k = 0; k < h; k++) {
                    T &x = a[j + k], &y = a[j + k + h];
                    if (type == OR) {
                        if (!inv) y += x; else y -= x;
                    } else if (type == AND) {
                        if (!inv) x += y; else x -= y;
                    } else { // XOR
                        T u = x, v = y;
                        x = u + v; y = u - v;
                        if (inv) { 
                            // 注意：如果是整数或浮点数直接 /2
                            // 如果是 ModInt，此处应乘上 2 的逆元
                            x /= 2; y /= 2; 
                        }
                    }
                }
            }
        }
    }

    // 卷积通用接口
    static vector<T> convolution(vector<T> a, vector<T> b, Type type) {
        int n = 1, sz = max(a.size(), b.size());
        while (n < sz) n <<= 1;
        a.resize(n, T(0)); b.resize(n, T(0));

        transform(a, type, false);
        transform(b, type, false);
        for (int i = 0; i < n; i++) a[i] *= b[i];
        transform(a, type, true);
        return a;
    }

    // 快捷调用
    static vector<T> convOR(vector<T> a, vector<T> b) { return convolution(a, b, OR); }
    static vector<T> convAND(vector<T> a, vector<T> b) { return convolution(a, b, AND); }
    static vector<T> convXOR(vector<T> a, vector<T> b) { return convolution(a, b, XOR); }
};
FWTSolver<ll> solver;

// 超集和
void zeta(vector<ll> &a) {
    for (int b = 0; b < M; b++) {
        int len = 1 << b;
        for (int i = 0; i < (1<<M); i += len << 1) {
            for (int j = 0; j < len; j++) {
                a[i + j] += a[i + j + len];
                if (a[i + j] >= MOD) a[i + j] -= MOD;
            }
        }
    }
}

// 其实就是zeta的逆变换
void mobius(vector<ll> &a) {
    for (int b = 0; b < M; b++) {
        int len = 1 << b;
        for (int i = 0; i < (1<<M); i += len << 1) {
            for (int j = 0; j < len; j++) {
                a[i + j] -= a[i + j + len];
                if (a[i + j] < 0) a[i + j] += MOD;
            }
        }
    }
}

void solve() {
    int n; cin >> n;
    vector<ll> cnt(1 << M, 0);
    for (int i = 0;i<n;i++) {
        ll x; cin >> x;
        cnt[x] ++;
    }

    vector<ll> fib(1 << M, 0);
    fib[1] = 1;
    for (int i = 2;i<(1<<M);i++) {
        fib[i] = (fib[i - 1] + fib[i - 2]) % MOD;
    }

    vector<ll> A(1 << M, 0);
    for (int i = 0;i<(1<<M);i++) {
        // 超酷枚举子集
        for (int j = i;; j = (j - 1) & i) {
            A[i] = (A[i] + cnt[j] * cnt[i ^ j]) % MOD;
            if (j == 0) break;
        }
        A[i] = A[i] * fib[i] % MOD;
    }

    vector<ll> C = solver.convXOR(cnt, cnt);
    for (int i = 0;i<(1<<M);i++) {
        C[i] = C[i] % MOD * fib[i] % MOD;
    }

    vector<ll> B = cnt;
    for (int i = 0;i<(1<<M);i++) {
        B[i] = B[i] * fib[i] % MOD;
    }

    zeta(A), zeta(B), zeta(C);
    vector<ll> T(1 << M, 0);
    for (int i = 0;i<(1<<M);i++) {
        T[i] = A[i] * B[i] % MOD * C[i] % MOD;
    }
    mobius(T);

    ll ans = 0;
    for (int i = 0;i<M;i++) {
        ans = (ans + T[1 << i]) % MOD;
    }

    cout << ans << endl;
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

// 这里还有一个十进制版本的 sosdp

// 数组长度必须为 10 的非负整数次幂。
void zeta(vector<int>& c, vector<int>& s, vector<int>& q) {
    const int V = c.size();
    assert(s.size() == c.size() && q.size() == c.size());
    int length = V;
    while (length > 1 && length % 10 == 0) length /= 10;
    assert(length == 1);
    for (int base = 1; base < V; base *= 10) {
        int step = base * 10;
        for (int l = 0; l < V; l += step) {
            for (int off = 0; off < base; off++) {
                for (int d = 8; d >= 0; d--) {
                    int i = l + d * base + off;
                    c[i] += c[i + base];

                    s[i] += s[i + base];
                    if (s[i] >= MOD) s[i] -= MOD;

                    q[i] += q[i + base];
                    if (q[i] >= MOD) q[i] -= MOD;
                }
            }
        }
    }
}

void mobius(vector<int>& f) {
    const int V = f.size();
    int length = V;
    while (length > 1 && length % 10 == 0) length /= 10;
    assert(length == 1);
    for (int base = 1; base < V; base *= 10) {
        int step = base * 10;
        for (int l = 0; l < V; l += step) {
            for (int off = 0; off < base; off++) {
                for (int d = 0; d < 9; d++) {
                    int i = l + d * base + off;
                    f[i] -= f[i + base];
                    if (f[i] < 0) f[i] += MOD;
                }
            }
        }
    }
}