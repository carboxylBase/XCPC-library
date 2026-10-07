#include <bits/stdc++.h>
using namespace std;
using ll = long long;

void solve() {
    string s, t;
    cin >> s >> t;
    char c = s[0];
    for (int i = 1; i < (int)s.size(); i++) {
        if (c == 'd') c = s[i];
        else if (c == s[i]) c = 'd';
        else c = 'a' + 'b' + 'c' - c - s[i];
    }
    vector<int> cnt(4, 0);
    ll ans = 0;
    for (char ch : t) {
        vector<int> ncnt(4, 0);
        int x = ch - 'a';
        ncnt[x] += cnt[3];
        for (int j = 0; j < 3; j++) {
            if (j == x) ncnt[3] += cnt[j];
            else ncnt[3 - j - x] += cnt[j];
        }
        ncnt[x]++;
        swap(cnt, ncnt);
        ans += cnt[c - 'a'];
    }
    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int t;
    cin >> t;
    while (t--) solve();
}
