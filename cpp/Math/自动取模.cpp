#include <iostream>

template <int Mod>
struct ModInt {
    int v;
    ModInt(long long _v = 0) { v = (_v % Mod + Mod) % Mod; }
    
    // 基础运算
    ModInt& operator+=(const ModInt& o) { v += o.v; if (v >= Mod) v -= Mod; return *this; }
    ModInt& operator-=(const ModInt& o) { v -= o.v; if (v < 0) v += Mod; return *this; }
    ModInt& operator*=(const ModInt& o) { v = 1LL * v * o.v % Mod; return *this; }
    
    // 快速幂与逆元
    ModInt pow(long long n) const {
        ModInt res(1), a(*this);
        while (n > 0) {
            if (n & 1) res *= a;
            a *= a; n >>= 1;
        }
        return res;
    }
    ModInt inv() const { return pow(Mod - 2); } // 费马小定理，要求 Mod 是质数
    
    ModInt& operator/=(const ModInt& o) { return *this *= o.inv(); }
    
    // 友元运算符
    friend ModInt operator+(ModInt a, const ModInt& b) { return a += b; }
    friend ModInt operator-(ModInt a, const ModInt& b) { return a -= b; }
    friend ModInt operator*(ModInt a, const ModInt& b) { return a *= b; }
    friend ModInt operator/(ModInt a, const ModInt& b) { return a /= b; }
    
    // 输入输出
    friend std::ostream& operator<<(std::ostream& os, const ModInt& a) { return os << a.v; }
    friend std::istream& operator>>(std::istream& is, ModInt& a) { long long v; is >> v; a = ModInt(v); return is; }
    
    bool operator==(const ModInt& o) const { return v == o.v; }
    bool operator!=(const ModInt& o) const { return v != o.v; }
};