#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll N = 2000000;

struct Dsu{
    int fa[N];
    void init(int n){
        for (int i = 1;i<n+1;i++){
            fa[i] = i;
        }
    }
    int find(int x){
        if (x == fa[x]){
            return x;
        }
        return fa[x] = find(fa[x]);
    }
    void merge(int x,int y){
        x = find(x),y = find(y);
        if (x == y){
            return;
        }
        fa[x] = y;
        return;
    }
};

struct dsu {
    struct edge {
        int x, y, sx;
        bool ok;
    };

    static const int N = 2000005;

    int fa[N], sz[N], top;
    edge st[N];

    void init(int n) {
        top = 0;
        for (int i = 1; i <= n; i++) {
            fa[i] = i;
            sz[i] = 1;
        }
    }

    int find(int x) {
        while (x != fa[x]) x = fa[x];
        return x;
    }

    bool merge(int u, int v) {
        int x = find(u), y = find(v);

        if (x == y) {
            st[++top] = {x, y, sz[x], false};
            return false;
        }

        if (sz[x] < sz[y]) swap(x, y);

        st[++top] = {x, y, sz[x], true};

        fa[y] = x;
        sz[x] += sz[y];

        return true;
    }

    void back() {
        edge e = st[top--];

        if (!e.ok) return;

        fa[e.y] = e.y;
        sz[e.x] = e.sx;
    }
};