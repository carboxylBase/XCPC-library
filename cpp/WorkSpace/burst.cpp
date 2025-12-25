#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <numeric>
#include <cstdio>
#include <cstdint>
#include <utility>
#include <functional>
using namespace std;
struct FS {
    static const int BS = 1 << 20;
    int i = 0, sz = 0;
    char b[BS];
    inline char rd() {
        if (i >= sz) {
            sz = (int)fread(b, 1, BS, stdin);
            i = 0;
            if (sz == 0) return 0;
        }
        return b[i++];
    }
    template<class T>
    bool rdInt(T &o) {
        char ch;
        do { ch = rd(); if (!ch) return false; } while (ch <= ' ');
        T s = 1;
        if (ch == '-') { s = -1; ch = rd(); }
        T x = 0;
        while (ch > ' ') {
            x = x * 10 + (ch - '0');
            ch = rd();
        }
        o = x * s;
        return true;
    }
    bool rdChar(char &o) {
        char ch;
        do { ch = rd(); if (!ch) return false; } while (ch <= ' ');
        o = ch;
        return true;
    }
};
inline static unsigned long long pP(int z, int v) {
    return ( (unsigned long long)z << 20 ) | (unsigned long long)v;
}
inline static unsigned long long pE(int z, int u, int v) {
    return ( (unsigned long long)z << 40 ) | ( (unsigned long long)u << 20 ) | (unsigned long long)v;
}
inline static void uE(unsigned long long k, int &z, int &u, int &v) {
    z = (int)(k >> 40);
    u = (int)((k >> 20) & ((unsigned long long)1 << 20) - 1);
    v = (int)(k & ((unsigned long long)1 << 20) - 1);
}
struct E {
    int x, y, z;
};
struct OR {
    unsigned char t;
    int a, b, sb;
    int z;
};
struct RDSU {
    int n;
    vector<int> p, s;
    vector<unsigned char> x;
    vector<OR> stk;
    vector<int> cf;
    int g;
    RDSU(int sz = 0, int Z = 0) { init(sz, Z); }
    void init(int sz, int Z) {
        n = sz;
        p.resize(n);
        s.assign(n, 1);
        x.assign(n, 0);
        iota(p.begin(), p.end(), 0);
        stk.clear();
        cf.assign(Z + 1, 0);
        g = Z;
    }
    inline pair<int,int> fd(int a) {
        int px = 0;
        while (p[a] != a) {
            px ^= x[a];
            a = p[a];
        }
        return {a, px};
    }
    inline void ut(int a, int b, int z) {
        auto [rx, px] = fd(a);
        auto [ry, py] = fd(b);
        if (rx == ry) {
            if ( (px ^ py) != 1 ) {
                int sv = ++cf[z];
                if (sv == 1) --g;
                stk.push_back({1, 0, 0, 0, z});
            }
            return;
        }
        if (s[rx] > s[ry]) {
            swap(rx, ry);
            swap(px, py);
        }
        stk.push_back({0, rx, ry, s[ry], 0});
        p[rx] = ry;
        x[rx] = (unsigned char)(px ^ py ^ 1);
        s[ry] += s[rx];
    }
    inline int sn() const { return (int)stk.size(); }
    inline void rb(int sp) {
        while ((int)stk.size() > sp) {
            auto op = stk.back();
            stk.pop_back();
            if (op.t == 0) {
                int rx = op.a, ry = op.b;
                p[rx] = rx;
                x[rx] = 0;
                s[ry] = op.sb;
            } else {
                int z = op.z;
                int sv = --cf[z];
                if (sv == 0) ++g;
            }
        }
    }
};
struct QO {
    char t;
    int u, v, z;
};
int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    FS fs;
    int n, m, Z, Q;
    if (!fs.rdInt(n)) return 0;
    fs.rdInt(m); fs.rdInt(Z); fs.rdInt(Q);
    vector<array<int,3>> ie;
    ie.reserve(m);
    int nt = 0;
    vector<QO> qs(Q + 1);
    vector<unsigned long long> pk;
    pk.reserve((size_t)2 * (m + Q) + 5);
    for (int j = 0; j < m; ++j) {
        int u, v, z;
        fs.rdInt(u); fs.rdInt(v); fs.rdInt(z);
        ie.push_back({u, v, z});
        pk.push_back(pP(z, u));
        pk.push_back(pP(z, v));
    }
    vector<char> iq(Q + 1, 0);
    for (int j = 1; j <= Q; ++j) {
        char ch;
        fs.rdChar(ch);
        qs[j].t = ch;
        if (ch == 'T') {
            int u, v, z;
            fs.rdInt(u); fs.rdInt(v); fs.rdInt(z);
            qs[j].u = u; qs[j].v = v; qs[j].z = z;
            pk.push_back(pP(z, u));
            pk.push_back(pP(z, v));
            ++nt;
        } else {
            iq[j] = 1;
        }
    }
    sort(pk.begin(), pk.end());
    pk.erase(unique(pk.begin(), pk.end()), pk.end());
    unordered_map<unsigned long long, int> pid;
    pid.reserve((size_t)(pk.size() * 1.3) + 10);
    pid.max_load_factor(0.7f);
    for (int j = 0; j < (int)pk.size(); ++j) pid[pk[j]] = j;
    int tp = (int)pk.size();
    vector<vector<E>> sg(4 * (Q + 5));
    function<void(int,int,int,int,int,const E&)> aS = [&](int id, int l, int r, int ql, int qr, const E &e) {
        if (ql <= l && r <= qr) {
            sg[id].push_back(e);
            return;
        }
        int mid = (l + r) >> 1;
        if (ql <= mid) aS(id<<1, l, mid, ql, qr, e);
        if (qr > mid) aS(id<<1|1, mid+1, r, ql, qr, e);
    };
    unordered_map<unsigned long long, int> stt;
    stt.reserve((size_t)((m + nt) * 1.3) + 10);
    stt.max_load_factor(0.7f);
    auto aI = [&](int l, int r, int u, int v, int z) {
        if (l > r) return;
        int ax = pid[pP(z, u)];
        int bx = pid[pP(z, v)];
        E e{ax, bx, z};
        aS(1, 1, Q, l, r, e);
    };
    for (auto &e : ie) {
        int u = e[0], v = e[1], z = e[2];
        unsigned long long k = pE(z, u, v);
        auto it = stt.find(k);
        if (it == stt.end()) stt.emplace(k, 1);
        else if (it->second == 0) it->second = 1;
    }
    for (int j = 1; j <= Q; ++j) {
        if (qs[j].t != 'T') continue;
        int u = qs[j].u, v = qs[j].v, z = qs[j].z;
        unsigned long long k = pE(z, u, v);
        auto it = stt.find(k);
        if (it == stt.end() || it->second == 0) {
            stt[k] = j;
        } else {
            int l = it->second;
            int r = j - 1;
            aI(l, r, u, v, z);
            it->second = 0;
        }
    }
    for (auto &kv : stt) {
        if (kv.second == 0) continue;
        int z, u, v;
        uE(kv.first, z, u, v);
        aI(kv.second, Q, u, v, z);
    }
    RDSU dsu(tp, Z);
    vector<int> ans;
    ans.reserve(Q);
    function<void(int,int,int)> d = [&](int id, int l, int r) {
        int sp = dsu.sn();
        for (auto &e : sg[id]) dsu.ut(e.x, e.y, e.z);
        if (l == r) {
            if (iq[l]) ans.push_back(dsu.g);
        } else {
            int mid = (l + r) >> 1;
            d(id<<1, l, mid);
            d(id<<1|1, mid+1, r);
        }
        dsu.rb(sp);
    };
    d(1, 1, Q);
    for (int x : ans) {
        cout << x << '\n';
    }
    return 0;
}