#include <bits/stdc++.h>
using namespace std;
using pii = pair<int,int>;

void solve() {
    int n,m,s; cin >> n >> m >> s;
    vector<int> a(n+1);
    for (int i=1; i<=n; i++) { cin >> a[i]; a[i]&=1; }
    vector<vector<pii>> g(n+1);
    for (int i=1; i<=m; i++) {
        int u,v; cin >> u >> v;
        g[u].push_back({v,i}); g[v].push_back({u,i});
    }
    vector<bool> usd(m+1),vis(n+1);
    auto dfs = [&](auto&& self,int u)->void {
        for (auto [v,idx]:g[u]) {
            if (vis[v]) continue;
            vis[v]=1; usd[idx]=1; self(self,v);
        }
    };
    vis[s]=1; dfs(dfs,s);
    for (int i=1; i<=n; i++) {
        if (!vis[i] && a[i]) { cout << "No\n"; return; }
    }
    int tot=0;
    for (int i=1; i<=n; i++) tot^=a[i];
    vector<int> huan;
    vector<bool> inhuan(n+1);
    for (int i=1; i<=n; i++) vis[i]=0;
    vector<int> stk,stk_pos(n+1,-1);
    auto dfs1 = [&](auto&& self,int u)->void {
        stk.push_back(u); stk_pos[u]=(int)stk.size()-1;
        if (!huan.empty()) return;
        for (auto [v,idx]:g[u]) {
            if (!huan.empty()) return;
            if (stk_pos[v]!=-1) {
                if (((stk_pos[u]-stk_pos[v])&1)==0) {
                    if (!huan.empty()) continue;
                    for (int i=stk_pos[v]; i<=stk_pos[u]; i++) huan.push_back(stk[i]);
                    return;
                }
            }
            if (!vis[v]) { vis[v]=1; self(self,v); }
        }
        stk_pos[u]=-1; stk.pop_back();
    };
    vis[s]=1; dfs1(dfs1,s);
    if ((tot&1) && huan.empty()) { cout << "No\n"; return; }
    for (auto v:huan) inhuan[v]=1;
    cout << "Yes\n";
    vector<int> ans;
    auto dfs2 = [&](auto&& self,int u,int fa)->void {
        if ((tot&1) && inhuan[u]) {
            tot=0;
            int pos=0;
            for (int i=0; i<(int)huan.size(); i++) if (huan[i]==u) pos=i;
            for (int i=0; i<(int)huan.size(); i++) {
                ans.push_back(huan[(pos+i+1)%huan.size()]);
                a[huan[(pos+i+1)%huan.size()]]^=1;
            }
        }
        for (auto [v,idx]:g[u]) {
            if (v==fa || usd[idx]==0) continue;
            ans.push_back(v); a[v]^=1;
            self(self,v,u);
            ans.push_back(u); a[u]^=1;
        }
        if ((a[u]&1) && u!=fa) {
            ans.push_back(fa); ans.push_back(u);
            a[u]^=1; a[fa]^=1;
        }
    };
    dfs2(dfs2,s,s);
    cout << ans.size() << '\n';
    for (auto v:ans) cout << v << ' ';
    cout << '\n';
}
int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int t; cin >> t;
    while (t--) solve();
}
