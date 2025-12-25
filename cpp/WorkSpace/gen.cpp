#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // freopen("input.txt", "w", stdout);

    int n = 3, m = 2, C = 2, Q = 5;
    cout << n << " " << m << " " << C << " " << Q << "\n";

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

    for (int i = 0; i < m; i++) {
        int u = rng() % n + 1;
        int v = rng() % n + 1;
        while (u == v) v = rng() % n + 1;
        if (u > v) swap(u, v);
        int c = rng() % C + 1;
        cout << u << " " << v << " " << c << "\n";
    }

    for (int i = 0; i < Q; i++) {
        if (rng() % 2) {
            cout << "Q\n";
        } else {
            int u = rng() % n + 1;
            int v = rng() % n + 1;
            while (u == v) v = rng() % n + 1;
            if (u > v) swap(u, v);
            int c = rng() % C + 1;
            cout << "T " << u << " " << v << " " << c << "\n";
        }
    }
    return 0;
}
