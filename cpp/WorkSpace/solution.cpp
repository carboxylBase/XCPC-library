#include <bits/stdc++.h>

using namespace std;
using ll = long long;

void solve() {
    int n; cin >> n;
    string s; cin >> s;
    vector<bool> usd(n, 0);
    vector<int> q1, q2;
    // 1. 消除同类型匹配
    for (int i = 0; i < n; i++) {
        if (s[i] == '(') {
            q1.push_back(i);
        } else if (s[i] == '[') {
            q2.push_back(i);
        } else if (s[i] == ')') {
            if (!q1.empty()) {
                usd[q1.back()] = 1;
                q1.pop_back();
                usd[i] = 1;
            }
        } else {
            if (!q2.empty()) {
                usd[q2.back()] = 1;
                q2.pop_back();
                usd[i] = 1;
            }
        }
    }

    string t;
    for (int i = 0; i < n; i++) {
        if (!usd[i]) t.push_back(s[i]);
    }

    int pairs_found = 0; // 记录在这个循环中找到了多少对
    int L = 0, R = 0;
    
    // 2. 消除 Open...Close 混合对
    for (auto v : t) {
        if (v == '(' || v == '[') {
            L++;
        } else {
            if (L > 0) {
                L--;
                pairs_found++; 
            } else {
                R++;
            }
        }
    }

    // 3. 计算最终代价
    int ans = pairs_found;
    ans += R / 2;
    ans += L / 2;

    if (R % 2 != 0) {
        // 如果剩余奇数个，说明有一对 Close...Open 没法配
        // 正常代价是 2
        // 但如果之前有 pairs_found，可以优化成 1
        if (pairs_found > 0) {
            ans += 1;
        } else {
            ans += 2;
        }
    }

    cout << ans << endl;
}

signed main() {freopen("input.txt", "r", stdin);
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int t = 1;
    cin >> t;

    while (t--) {
        solve();
    }

    return 0;
}