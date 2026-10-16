#include <bits/stdc++.h>
using namespace std;
using ll = long long;

// 默认模数可通过第二个模板参数修改：Mat<ll, 998244353>。
// 元素使用整数类型，保持在 [0, Mod)；Mod > 0。
template<typename T, ll Mod = 1000000007>
struct Mat {
    static_assert(Mod > 0);
    int n, m;
    vector<vector<T>> a;
    Mat(int n = 0, int m = 0) : n(n), m(m), a(n, vector<T>(m, 0)) {}
    static Mat identity(int n) {
        Mat ret(n, n);
        for (int i = 0; i < n; ++i) ret.a[i][i] = 1 % Mod;
        return ret;
    }
    Mat operator+(const Mat& b) const {
        assert(n == b.n && m == b.m);
        Mat ret(n, m);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                ret.a[i][j] = ((__int128)a[i][j] + b.a[i][j]) % Mod;
        return ret;
    }
    Mat operator*(const Mat& b) const {
        assert(m == b.n);
        Mat ret(n, b.m);
        for (int i = 0; i < n; ++i)
            for (int k = 0; k < m; ++k)
                for (int j = 0; j < b.m; ++j)
                    ret.a[i][j] = (ret.a[i][j] + (__int128)a[i][k] * b.a[k][j]) % Mod;
        return ret;
    }
    Mat& operator+=(const Mat& b) { return *this = *this + b; }
    Mat& operator*=(const Mat& b) { return *this = *this * b; }
};
template<typename T, ll Mod>
Mat<T, Mod> qpow(Mat<T, Mod> a, ll b) {
    assert(a.n == a.m && b >= 0);
    auto ret = Mat<T, Mod>::identity(a.n);
    for (; b; b >>= 1, a *= a)
        if (b & 1) ret *= a;
    return ret;
}
