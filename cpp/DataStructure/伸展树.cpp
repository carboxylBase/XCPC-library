// It's a wonderful life.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define DEBUG 1
const ll N = 2000000;
const ll MOD = 998244353;
const ll MAX = 1e18;

/*
可重集合
rank 得到的是比 val 小的 Node 的数量
*/
template<typename T, typename Compare = std::less<T>>
struct Splay{
    int rt,tot;
    struct Node{
        T val;
        int cnt,fa,ch[2],siz;
    }nodes[N];

    Compare cmp;
    bool lt(const T &a, const T &b) const { return cmp(a,b); }
    bool eq(const T &a, const T &b) const { return !cmp(a,b) && !cmp(b,a); }

    int newNode(const T &x){
        ++tot;
        nodes[tot].val = x;
        nodes[tot].cnt = nodes[tot].siz = 1;
        nodes[tot].fa = nodes[tot].ch[0] = nodes[tot].ch[1] = 0;
        return tot;
    }

    void init(){
        // 初始化 nodes[0]
        nodes[0].cnt = nodes[0].siz = nodes[0].fa = nodes[0].ch[0] = nodes[0].ch[1] = 0;
        tot = 0;
        int a = newNode(T(-INT_MAX));
        int b = newNode(T(INT_MAX));
        rt = a;
        nodes[rt].ch[1] = b;
        nodes[b].fa = rt;
        update(rt);
    }

    void update(int x){
        nodes[x].siz = nodes[nodes[x].ch[0]].siz + nodes[nodes[x].ch[1]].siz + nodes[x].cnt;
    }

    void rot_left(int x){
        int y = nodes[x].fa,z = nodes[y].fa;
        nodes[y].ch[1] = nodes[x].ch[0]; nodes[nodes[x].ch[0]].fa = y;
        nodes[x].ch[0] = y; nodes[y].fa = x;
        nodes[z].ch[nodes[z].ch[1] == y] = x; nodes[x].fa = z;
        update(y); update(x);
    }

    void rot_right(int x){
        int y = nodes[x].fa,z = nodes[y].fa;
        nodes[y].ch[0] = nodes[x].ch[1]; nodes[nodes[x].ch[1]].fa = y;
        nodes[x].ch[1] = y; nodes[y].fa = x;
        nodes[z].ch[nodes[z].ch[1] == y] = x; nodes[x].fa = z;
        update(y); update(x);
    }

    int getlr(int x){
        return nodes[nodes[x].fa].ch[1] == x;
    }

    void rotate(int x){
        if (getlr(x)) rot_left(x); else rot_right(x);
    }

    void splay(int x,int target){
        if (!target) rt = x;
        while (nodes[x].fa != target){
            int y = nodes[x].fa, z = nodes[y].fa;
            if (z != target){
                if (getlr(x) == getlr(y)) rotate(y);
                else rotate(x);
            }
            rotate(x);
        }
    }

    void find(const T &x){
        if (!rt) return;
        int p = rt;
        while (!eq(nodes[p].val, x) && nodes[p].ch[ lt(nodes[p].val, x) ? 1 : 0 ]){
            p = nodes[p].ch[ lt(nodes[p].val, x) ? 1 : 0 ];
        }
        splay(p,0);
    }

    int pre(const T &x){
        find(x);
        if (lt(nodes[rt].val, x)) return rt;
        int p = nodes[rt].ch[0];
        while (nodes[p].ch[1]) p = nodes[p].ch[1];
        splay(p,0);
        return p;
    }

    int suc(const T &x){
        find(x);
        if (lt(x, nodes[rt].val)) return rt;
        int p = nodes[rt].ch[1];
        while (nodes[p].ch[0]) p = nodes[p].ch[0];
        splay(p,0);
        return p;
    }

    void insert(const T &x){
        int p = rt, fp = 0;
        while (p && !eq(nodes[p].val, x)){
            fp = p; p = nodes[p].ch[ lt(nodes[p].val, x) ? 1 : 0 ];
        }
        if (!p){
            p = newNode(x);
            nodes[fp].ch[ lt(nodes[fp].val, x) ? 1 : 0 ] = p;
            nodes[p].fa = fp;
        }else{
            nodes[p].cnt++;
        }
        splay(p,0);
    }

    void del(const T &x){
        int xPre = pre(x), xSuc = suc(x);
        splay(xPre,0); splay(xSuc,xPre);
        int d = nodes[xSuc].ch[0];
        if (--nodes[d].cnt){
            splay(d,0);
        }else{
            nodes[xSuc].ch[0] = 0;
            update(xSuc); update(xPre);
        }
    }

    int kth(long long x){
        int p = rt;
        while (1){
            int v = nodes[p].ch[0];
            if (nodes[v].siz + nodes[p].cnt < x){
                x -= nodes[v].siz + nodes[p].cnt;
                p = nodes[p].ch[1];
            }else{
                if (nodes[v].siz < x){
                    splay(p,0);
                    return p;
                }else p = v;
            }
        }
    }

    int getRank(const T &x){
        find(x);
        return nodes[nodes[rt].ch[0]].siz;
    }
};

struct Node {
    ll val;
    Node(ll VAL_ = 0) {val = VAL_;}
    bool operator<(const Node& A) const {
        return val < A.val;
    }
};
Splay<Node> solver;

void solve(){
    int n,m;cin >> n >> m;
    solver.init();
    for (int i = 1;i<n+1;i++){
        int x;cin >> x;solver.insert(x);
    }
    int last = 0,ans = 0;
    while (m--){
        int opt,x;cin >> opt >> x;
        x ^= last;
        if (opt == 1){
            solver.insert(x);
        }else if (opt == 2){
            solver.del(x);
        }else if (opt == 3){
            solver.insert(x);
            last = solver.getRank(x);
            solver.del(x);
        }else if (opt == 4){
            last = solver.nodes[solver.kth(x+1)].val.val;
        }else if (opt == 5){
            last = solver.nodes[solver.pre(x)].val.val;
        }else{
            last = solver.nodes[solver.suc(x)].val.val;
        }
        if (opt > 2){
            ans ^= last;
        }
    }

    cout << ans;
    return;
}

signed main(){
    freopen("input.txt","r",stdin);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int _ = 1;
    // cin >> _;
    while (_--){
        solve();
    }
    return 0;
}