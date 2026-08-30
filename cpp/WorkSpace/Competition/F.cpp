// It's a wonderful life.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
#define DEBUG 1
const ll N = 2000000;
const ll MOD = 998244353;
const ll INF = 1e18;

int n,a[N],b[N];
vector<pair<int,int>> ans;

void solve() {
    cin >> n;
    for (int i = 1;i<n+1;i++) {
        cin >> a[i];
        b[i] = i;
    }

    for (int i = 1;i<n+1;i++) {
        int cur;
        for (int j = i;j<n+1;j++) {
            if (b[j] == a[i]) {
                cur = j;break;
            }
        }

        while (cur + 1 <= n) {
            // cout << b[cur + 1] << " " << b[cur] << endl;
            ans.push_back({b[cur + 1],b[cur]});
            swap(b[cur],b[cur + 1]);
            cur++;
        }

        cur = n;
        while (cur - 1 >= i) {
            // cout << b[cur] << " " << b[cur - 1] << endl;
            ans.push_back({b[cur],b[cur - 1]});
            swap(b[cur],b[cur - 1]);
            cur--;
        }
    }
    cout << ans.size() << endl;
    for (auto v : ans) {
        cout << v.first << " " << v.second << endl;
    }
    return; 
}

signed main() {
    //freopen("input.txt","r",stdin);
    ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
    int _ = 1;
    // cin >> _;
    while (_--){
        solve();
    }
    return 0;
}
