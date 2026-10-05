#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using pii = pair<int, int>;

struct Primes {
    bool notPrime[100001]{};
    vector<int> primes;
    void sieve(int maxn) {
        for (ll i = 2; i <= maxn; i++) {
            if (!notPrime[i]) primes.push_back(i);
            for (int p : primes) {
                if (i * p > maxn) break;
                notPrime[i * p] = true;
                if (i % p == 0) break;
            }
        }
    }
} solver;

void solve() {
    ll n;
    cin >> n;
    ll m = n;
    vector<pii> a;
    for (int v : solver.primes) {
        if (n % v == 0) {
            pii z(v, 0);
            while (n % v == 0) {
                n /= v;
                z.second++;
            }
            a.push_back(z);
        }
    }
    if (n > 1) a.emplace_back((int)n, 1);

    int phi_n = m;
    for (auto [v, w] : a) {
        phi_n = 1LL * phi_n * (v - 1) / v;
    }

    vector<int> divisors;
    for (int d = 1; 1LL * d * d <= phi_n; d++) {
        if (phi_n % d != 0) continue;
        divisors.push_back(d);
        if (d != phi_n / d) divisors.push_back(phi_n / d);
    }
    sort(divisors.begin(), divisors.end());

    unordered_map<int, ll> dp;
    for (int d : divisors) dp[d] = 1;
    for (auto [p, e] : a) {
        if (p == 2 && e >= 3) {
            ll cycle = 1LL << (e - 2);
            for (int d : divisors) {
                dp[d] *= gcd((ll)d, 2LL) * gcd((ll)d, cycle);
            }
        } else {
            ll q = 1;
            for (int i = 0; i < e; i++) q *= p;
            ll phi_q = q / p * (p - 1);
            for (int d : divisors) dp[d] *= gcd((ll)d, phi_q);
        }
    }

    vector<int> b;
    for (int v : solver.primes) {
        if (phi_n % v == 0) {
            while (phi_n % v == 0) phi_n /= v;
            b.push_back(v);
        }
    }
    if (phi_n > 1) b.push_back(phi_n);

    reverse(divisors.begin(), divisors.end());
    for (int p : b) {
        for (int d : divisors) {
            if (d % p == 0) dp[d] -= dp[d / p];
        }
    }

    ll ans = 0;
    for (auto [v, w] : dp) ans += w * v;
    cout << ans << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solver.sieve(100000);
    int t;
    cin >> t;
    while (t--) solve();
}
