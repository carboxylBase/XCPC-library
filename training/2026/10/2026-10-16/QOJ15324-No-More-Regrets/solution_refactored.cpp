#include <bits/stdc++.h>
using namespace std;
using i64 = long long;
using u64 = unsigned long long;
constexpr int B = 1500;

struct Stats {
    u64 mn = 0, mx = 0, prod = 0;
    Stats operator-(const Stats& b) const {
        return {mn - b.mn, mx - b.mx, prod - b.prod};
    }
};

struct Block {
    vector<i64> a, mn, mx;
    vector<Stats> pre;
    i64 same = -1, add = 0;

    explicit Block(vector<i64> values) : a(move(values)),
        mn(a.size()), mx(a.size()), pre(a.size() + 1) { rebuild(); }

    int size() const { return (int)a.size(); }
    i64 get(int i) const { return (same == -1 ? a[i] : same) + add; }

    void rebuild() {
        pre[0] = {};
        for (int i = 0; i < size(); ++i) {
            mn[i] = i ? min(mn[i - 1], a[i]) : a[i];
            mx[i] = i ? max(mx[i - 1], a[i]) : a[i];
            pre[i + 1] = {pre[i].mn + u64(mn[i]),
                          pre[i].mx + u64(mx[i]),
                          pre[i].prod + u64(mn[i]) * u64(mx[i])};
        }
    }

    // All local intervals are half-open: [l, r).
    void modify(int l, int r, int op, i64 v) {
        if (l == 0 && r == size()) {
            if (op == 1) add += v;
            else same = v, add = 0;
            return;
        }
        if (same != -1 || add != 0) {
            for (int i = 0; i < size(); ++i) a[i] = get(i);
            same = -1; add = 0;
        }
        for (int i = l; i < r; ++i)
            a[i] = op == 1 ? a[i] + v : v;
        rebuild();
    }

    // Statistics for the first k prefix extrema, including lazy tags.
    Stats prefix(int k) const {
        if (same != -1) {
            u64 v = u64(same + add);
            return {v * k, v * k, v * v * k};
        }
        auto s = pre[k];
        u64 d = u64(add); // Arithmetic below is modulo 2^64, even for negative add.
        return {s.mn + d * k, s.mx + d * k,
                s.prod + d * (s.mn + s.mx) + d * d * k};
    }

    // Append this entire block to a query with current extrema lo, hi.
    u64 consume(i64& lo, i64& hi) const {
        if (same != -1) {
            i64 v = same + add;
            lo = min(lo, v); hi = max(hi, v);
            return u64(lo) * u64(hi) * size();
        }
        int p = partition_point(mn.begin(), mn.end(),
            [&](i64 x) { return x + add >= lo; }) - mn.begin();
        int q = partition_point(mx.begin(), mx.end(),
            [&](i64 x) { return x + add <= hi; }) - mx.begin();

        int l = min(p, q), r = max(p, q);
        auto left = prefix(l), right = prefix(r);
        auto middle = right - left;
        u64 ans = u64(lo) * u64(hi) * l;
        ans += p < q ? middle.mn * u64(hi) : middle.mx * u64(lo);
        ans += (prefix(size()) - right).prod;
        lo = min(lo, mn.back() + add);
        hi = max(hi, mx.back() + add);
        return ans;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n, queries;
    cin >> n >> queries;
    vector<i64> a(n);
    for (auto& x : a) cin >> x;
    vector<Block> blocks;
    for (int l = 0; l < n; l += B)
        blocks.emplace_back(vector<i64>(a.begin() + l, a.begin() + min(n, l + B)));

    while (queries--) {
        int op, l, r;
        cin >> op >> l >> r;
        --l; // Input [l, r] becomes zero-based [l-1, r).
        i64 v = 0;
        if (op != 3) cin >> v;
        i64 lo = LLONG_MAX, hi = 0;
        u64 ans = 0;
        for (int id = l / B; id <= (r - 1) / B; ++id) {
            auto& b = blocks[id];
            int L = max(l - id * B, 0), R = min(r - id * B, b.size());
            if (op != 3) b.modify(L, R, op, v);
            else if (L == 0 && R == b.size()) ans += b.consume(lo, hi);
            else for (int i = L; i < R; ++i) {
                i64 x = b.get(i);
                lo = min(lo, x); hi = max(hi, x);
                ans += u64(lo) * u64(hi);
            }
        }
        if (op == 3) cout << ans << '\n';
    }
}
