#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using pll = pair<ll,ll>;
vector<vector<pll>> g;
void solve() {
    int n,m,k; cin >> n >> m >> k;
    g.clear(); g.resize(n+1);
    vector<array<ll,3>> e(m);
    for (int i = 0; i < m; i++) {
        int x,y,z; cin >> x >> y >> z;
        g[x].push_back({y,z});
        e[i] = {z,x,y};
    }
    for (int i = 1; i <= n; i++) {
        sort(g[i].begin(), g[i].end(), [](pll l, pll r) { return l.second < r.second; });
    }
    ll tot = 0;
    vector<ll> ans;
    sort(e.begin(), e.end());
    struct Node { int nxt,pos; ll tot; };
    auto dfs = [&](auto&& self, vector<Node>& b) -> void {
        if (tot == k) return;
        vector<Node> c;
        for (auto v : b) {
            tot++;
            ans.push_back(v.tot);
            if (tot == k) return;
        }
        struct Z {
            int pos,cur;
            bool operator<(const Z& z) const {
                return g[pos][cur].second > g[z.pos][z.cur].second;
            }
        };
        priority_queue<Z> zs;
        for (auto v : b) {
            Z z{v.pos,0};
            if (z.cur < (int)g[z.pos].size()) zs.push(z);
        }
        while (!zs.empty() && tot + (ll)c.size() < k) {
            Z z = zs.top(); zs.pop();
            Node x;
            x.nxt = g[z.pos][z.cur].second;
            x.pos = g[z.pos][z.cur].first;
            x.tot = b[0].tot + 1;
            if (!c.empty() && x.nxt != c.back().nxt) {
                self(self,c);
                c.clear();
                if (tot == k) return;
                c.push_back(x);
            } else c.push_back(x);
            z.cur++;
            if (z.cur < (int)g[z.pos].size()) zs.push(z);
        }
        if (!c.empty()) self(self,c);
    };
    for (int i = 0; i < (int)e.size(); i++) {
        int j = i;
        for (; j < (int)e.size(); j++) if (e[j][0] != e[i][0]) break;
        j--;
        vector<Node> b;
        for (int p = i; p <= j; p++) b.push_back({(int)e[p][0],(int)e[p][2],1});
        dfs(dfs,b);
        i = j;
    }
    while ((int)ans.size() < k) ans.push_back(-1);
    for (auto v : ans) cout << v << '\n';
}
int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    solve();
}
