// It's a wonderful life.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define DEBUG 1
const ll N = 2000000;
const ll MOD = 998244353;
const ll INF = 1e18;

int n,m;
struct Node {
    int in;
    string opt,res;
    vector<int> son;
} nodes[N];

void dfs(int sta,deque<int> t[2],deque<int> f[2],vector<int> ans) {
    deque<int> nt[2],nf[2];
    vector<int> tag;

    while (!f[sta].empty()) {
        int u = f[sta][0];f[sta].pop_front();
        ans.push_back(u);
        for (auto v : nodes[u].son) {
            nodes[v].in--;
            tag.push_back(v);
            if (nodes[v].in == 0) {
                if (nodes[v].res == "true") {
                    if (nodes[v].opt == "set") {
                        t[1].push_back(v);
                    } else {
                        t[0].push_back(v);
                    }
                } else {
                    if (nodes[v].opt == "set") {
                        f[1].push_back(v);
                    } else {
                        f[0].push_back(v);
                    }
                }
            }
        }
    }

    // for (auto v : ans) {
    //     cout << v << " ";
    // }
    // cout << endl;

    for (int i = 0;i<t[1^sta].size();i++) {
        vector<int> n_ans;
        n_ans = ans;
        n_ans.push_back(t[1^sta][i]);
        deque<int> n_f[2],n_t[2];
        n_f[0] = f[0],n_f[1] = f[1];
        n_t[sta] = t[sta];
        for (int j = 0;j<t[1^sta].size();j++) {
            if (j!=i) n_t[1^sta].push_back(t[1^sta][j]);
        }

        for (auto v : nodes[t[1^sta][i]].son) {
            nodes[v].in--;
            if (nodes[v].in == 0) {
                if (nodes[v].res == "true") {
                    if (nodes[v].opt == "set") {
                        n_t[1].push_back(v);
                    } else {
                        n_t[0].push_back(v);
                    }
                } else {
                    if (nodes[v].opt == "set") {
                        n_f[1].push_back(v);
                    } else {
                        n_f[0].push_back(v);
                    }
                }
            }
        }

        dfs(sta^1,n_t,n_f,n_ans);

        for (auto v : nodes[t[1^sta][i]].son) {
            nodes[v].in++;
        }
    }

    for (auto v : tag) {
        nodes[v].in++;
    }

    if (ans.size() == n) {
        for (auto v : ans) {
            cout << v << " ";
        }
        cout << endl;
        exit(0);
    }
    return;
}

void solve() {
    cin >> n;
    for (int i = 1;i<n+1;i++) {
        cin >> nodes[i].opt >> nodes[i].res;
    }

    cin >> m;
    for (int i = 1;i<m+1;i++) {
        int u,v;cin >> u >> v;
        nodes[u].son.push_back(v);
        nodes[v].son.push_back(u);
        nodes[v].in++;
    }

    vector<int> ans;
    deque<int> f[2],t[2];
    for (int i = 1;i<n+1;i++) {
        if (!nodes[i].in) {
            if (nodes[i].res == "true") {
                if (nodes[i].opt == "set") {
                    t[1].push_back(i);
                } else {
                    t[0].push_back(i);
                }
            } else {
                if (nodes[i].opt == "set") {
                    f[1].push_back(i);
                } else {
                    f[0].push_back(i);
                }
            }
        }
    }

    // cout << f[0][0] << endl;

    dfs(0,t,f,ans);

    cout << -1 << endl;
    return;
}

signed main() {
    // freopen("input.txt","r",stdin);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int _ = 1;
    // cin >> _;
    while (_--){
        solve();
    }
    return 0;
}
