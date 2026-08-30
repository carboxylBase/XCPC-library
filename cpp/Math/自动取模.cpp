#include <iostream>
using ll = long long;

template<int P>
struct MODint {
    int v, c;
    constexpr MODint() : v(0), c(0) {}
    constexpr MODint(ll x) : v(0), c(0) {
        if (x != 0) {
            while (x % P == 0) x /= P, c++;
            v = x % P;
            if (v < 0) v += P;
        }
    }
    constexpr MODint(int v, int c) : v(v), c(c) {}
    constexpr int val() const {
        return c > 0 ? 0 : v;
    }
    constexpr int inv(int a) const {
        int b = P - 2, res = 1;
        for (; b; b >>= 1, a = 1LL * a * a % P) {
            if (b & 1) res = 1LL * res * a % P;
        }
        return res;
    }
    constexpr MODint& operator*=(MODint rhs) {
        if (v == 0 || rhs.v == 0) return *this = MODint();
        v = 1LL * v * rhs.v % P;
        c += rhs.c;
        return *this;
    }
    constexpr MODint& operator/=(MODint rhs) {
        v = 1LL * v * inv(rhs.v) % P;
        c -= rhs.c;
        return *this;
    }
    friend constexpr MODint operator*(MODint lhs, MODint rhs) {
        MODint res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MODint operator/(MODint lhs, MODint rhs) {
        MODint res = lhs;
        res /= rhs;
        return res;
    }
    friend std::istream& operator>>(std::istream& is, MODint& a) {
        ll x;
        is >> x;
        a = MODint(x);
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, const MODint& a) {
        return os << a.val();
    }
};