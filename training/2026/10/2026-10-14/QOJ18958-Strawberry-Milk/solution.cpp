#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 998244353;

namespace poly {

    ll qpow(ll base, ll k) {
        if (base == 0) return 0;
        ll res = 1;
        base %= MOD; base = (base + MOD) % MOD;
        k %= (MOD - 1); k = (k + MOD - 1) % (MOD - 1);
        while (k) {
            if (k & 1) {
                res = res * base % MOD;
            }
            base = base * base % MOD;
            k >>= 1;
        }
        return res;
    }

    ll inv(ll x) {
        return qpow(x, MOD - 2);
    }

    mt19937 MT;

    using Poly = vector<ll>;
    #define lg(x) ((x) == 0 ? -1 : __lg(x))
    #define Size(x) int(x.size())

    namespace NTT_ns {
        const long long G = 3, invG = inv(G);
        vector<int> rev;

        void NTT(ll* F, int len, int sgn) {
            rev.resize(len);
            for (int i = 1; i < len; ++i) {
                rev[i] = (rev[i >> 1] >> 1) | ((i & 1) * (len >> 1));
                if (i < rev[i]) swap(F[i], F[rev[i]]);
            }

            for (int tmp = 1; tmp < len; tmp <<= 1) {
                ll w1 = qpow(sgn ? G : invG, (MOD - 1) / (tmp << 1));
                for (int i = 0; i < len; i += tmp << 1) {
                    for (ll j = 0, w = 1; j < tmp; ++j, w = w * w1 % MOD) {
                        ll x = F[i + j];
                        ll y = F[i + j + tmp] * w % MOD;
                        F[i + j]         = (x + y) % MOD;
                        F[i + j + tmp]   = (x - y + MOD) % MOD;
                    }
                }
            }

            if (sgn == 0) {
                ll inv_len = inv(len);
                for (int i = 0; i < len; ++i)
                    F[i] = F[i] * inv_len % MOD;
            }
        }
    }

    using NTT_ns::NTT;

    Poly operator*(Poly F, Poly G) {
        int siz = Size(F) + Size(G) - 1;
        int len = 1 << (lg(siz - 1) + 1);

        if (siz <= 300) {
            Poly H(siz);
            for (int i = Size(F) - 1; ~i; --i)
                for (int j = Size(G) - 1; ~j; --j)
                    H[i + j] = (H[i + j] + F[i] * G[j]) % MOD;
            return H;
        }

        F.resize(len);
        G.resize(len);

        NTT(F.data(), len, 1);
        NTT(G.data(), len, 1);

        for (int i = 0; i < len; ++i)
            F[i] = F[i] * G[i] % MOD;

        NTT(F.data(), len, 0);
        F.resize(siz);
        return F;
    }

    Poly operator+(Poly F, Poly G) {
        int siz = max(Size(F), Size(G));
        F.resize(siz);
        G.resize(siz);
        for (int i = 0; i < siz; ++i)
            F[i] = (F[i] + G[i]) % MOD;
        return F;
    }

    Poly operator-(Poly F, Poly G) {
        int siz = max(Size(F), Size(G));
        F.resize(siz);
        G.resize(siz);
        for (int i = 0; i < siz; ++i)
            F[i] = (F[i] - G[i] + MOD) % MOD;
        return F;
    }

    Poly lsh(Poly F, int k) {
        F.resize(Size(F) + k);
        for (int i = Size(F) - 1; i >= k; --i)
            F[i] = F[i - k];
        for (int i = 0; i < k; ++i)
            F[i] = 0;
        return F;
    }

    Poly rsh(Poly F, int k) {
        int siz = Size(F) - k;
        for (int i = 0; i < siz; ++i)
            F[i] = F[i + k];
        F.resize(siz);
        return F;
    }

    Poly cut(Poly F, int len) {
        F.resize(len);
        return F;
    }

    Poly der(Poly F) {
        int siz = Size(F) - 1;
        for (int i = 0; i < siz; ++i)
            F[i] = F[i + 1] * (i + 1) % MOD;
        F.pop_back();
        return F;
    }

    Poly inte(Poly F) {
        F.emplace_back(0);
        for (int i = Size(F) - 1; i; --i)
            F[i] = F[i - 1] * inv(i) % MOD;
        F[0] = 0;
        return F;
    }

    Poly inv(Poly F) {
        int siz = Size(F);
        Poly G{inv(F[0])};
        for (int i = 2; (i >> 1) < siz; i <<= 1) {
            G = G + G - G * G * cut(F, i);
            G.resize(i);
        }
        G.resize(siz);
        return G;
    }

    Poly ln(Poly F) {
        return cut(inte(cut(der(F) * inv(F), Size(F))), Size(F));
    }

    Poly exp(Poly F) {
        int siz = Size(F);
        Poly G{1};
        for (int i = 2; (i >> 1) < siz; i <<= 1) {
            G = G * (Poly{1} - ln(cut(G, i)) + cut(F, i));
            G.resize(i);
        }
        G.resize(siz);
        return G;
    }
};

using namespace poly;

void solve() {
    int n, k; cin >> n >> k;
    
    if (n == 1) {
        cout << 0 << endl;
        return;
    }

    auto cal = [&](auto&&self, int l, int r) -> Poly {
        if (l == r) {
            Poly a(3);
            a[0] = 1, a[1] = l - 1, a[2] = 1;
            return a;
        }
        int mid = l + r >> 1;
        Poly L = self(self, l, mid), R = self(self, mid + 1, r);
        return L * R;
    };
    Poly res = cal(cal, 1, n - 1);
    ll ans = 0;
    for (int i = 0;i<res.size();i++) {
        if (i - n + 1 >= k) {
            ans = (ans + res[i]) % MOD;
        }
    }
    cout << ans << endl;
    return;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    cin >> t;
    while (t--) solve();
    return 0;
}
