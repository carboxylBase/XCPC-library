#include <bits/stdc++.h>
using namespace std;
using ll = unsigned long long;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using db = long double;
using pdd = pair<db, db>;
using i128 = __int128_t;

const ll N = 2000000;
const ll INF = 5e18;
const ll MOD = 1e9 + 7;
const ll B = 1500;

void solve() {
    int n, q; cin >> n >> q;
    vector<ll> a(n, 0);
    for (auto &v : a) cin >> v;
    while (n % B != 0) {
        n ++;
        a.push_back(0);
    }

    struct Block {
        long long laz_same;
        ll laz_add;
        ll mx[B], mn[B], pre_mn[B], pre_mx[B], pre_val[B], ord[B];
    } bs[n / B];

    auto upd = [&](int i) -> void {
        if (bs[i].laz_same != -1) {
            for (int j = 0;j<B;j++) {
                bs[i].ord[j] = bs[i].laz_same;
            }
            bs[i].laz_same = -1;
        }
        if (bs[i].laz_add != 0) {
            for (int j = 0;j<B;j++) {
                bs[i].ord[j] += bs[i].laz_add;
            }
            bs[i].laz_add = 0;
        }
        for (int j = 0;j<B;j++) {
            bs[i].mx[j] = bs[i].mn[j] = bs[i].ord[j];
            if (j > 0) {
                bs[i].mx[j] = max(bs[i].mx[j], bs[i].mx[j - 1]);
                bs[i].mn[j] = min(bs[i].mn[j], bs[i].mn[j - 1]);
            }
            bs[i].pre_mn[j] = bs[i].mn[j];
            bs[i].pre_mx[j] = bs[i].mx[j];
            bs[i].pre_val[j] = bs[i].mn[j] * bs[i].mx[j];
            if (j > 0) {
                bs[i].pre_mn[j] += bs[i].pre_mn[j - 1];
                bs[i].pre_mx[j] += bs[i].pre_mx[j - 1];
                bs[i].pre_val[j] += bs[i].pre_val[j - 1];
            }
        }
        return;
    };

    for (int i = 0;i<n/B;i++) {
        bs[i].laz_same = -1;
        bs[i].laz_add = 0;
        for (int j = 0;j<B;j++) {
            bs[i].ord[j] = a[i * B + j];
        }
        upd(i);
    }

    while (q --) {
        int opt, l, r; cin >> opt >> l >> r;
        l --, r --;
        if (opt == 1) {
            ll v; cin >> v;
            upd(l / B);
            if (r / B != l / B) {
                upd(r / B);
            }
            if (l / B + 1 >= r / B) {
                for (int i = l;i<=r;i++) {
                    bs[i / B].ord[i % B] += v;
                }
                upd(l / B);
                if (l / B != r / B) upd(r / B);
            } else {
                int z = l / B;
                for (;l%B!=0;l++) {
                    bs[l / B].ord[l % B] += v;
                }
                upd(z);

                for (;l/B!=r/B;l+=B) {
                    bs[l / B].laz_add += v;
                }
                for (int i = l;i<=r;i++) {
                    bs[l / B].ord[i % B] += v;
                }
                upd(l / B);
            }
        } else if (opt == 2) {
            ll v; cin >> v;
            upd(l / B);
            if (r / B != l / B) {
                upd(r / B);
            }
            if (l / B + 1 >= r / B) {
                for (int i = l;i<=r;i++) {
                    bs[i / B].ord[i % B] = v;
                }
                upd(l / B);
                if (l / B != r / B) upd(r / B);
            } else {
                int z = l / B;
                for (;l%B!=0;l++) {
                    bs[l / B].ord[l % B] = v;
                }
                upd(z);

                for (;l/B!=r/B;l+=B) {
                    bs[l / B].laz_same = v;
                    bs[l / B].laz_add = 0;
                }

                for (int i = l;i<=r;i++) {
                    bs[l / B].ord[i % B] = v;
                }
                upd(l / B);
            }
        } else {
            ll ans = 0;
            if (l / B + 1 >= r / B) {
                upd(l / B);
                if (r / B != l / B) {
                    upd(r / B);
                }
                ll mn = INF, mx = 0;
                for (int i = l;i<=r;i++) {
                    mn = min(mn, bs[i / B].ord[i % B]);
                    mx = max(mx, bs[i / B].ord[i % B]);
                    ans += mn * mx;
                }
            } else {
                ll mn = INF, mx = 0;
                upd(l / B);
                upd(r / B);
                while (l % B != 0) {
                    mn = min(mn, bs[l / B].ord[l % B]);
                    mx = max(mx, bs[l / B].ord[l % B]);
                    ans += mn * mx;
                    l ++;
                }
                while (l / B != r / B) {
                    int mn_pos = -1, mx_pos = -1;
                    int lo = 0, hi = B - 1;
                    while (lo <= hi) {
                        int mid = lo + hi >> 1;
                        ll res = bs[l / B].mn[mid];
                        if (bs[l / B].laz_same != -1) {
                            res = bs[l / B].laz_same;
                        }
                        if (bs[l / B].laz_add != 0) {
                            res += bs[l / B].laz_add;
                        }
                        if (res >= mn) {
                            mn_pos = mid;
                            lo = mid + 1;
                        } else {
                            hi = mid - 1;
                        }
                    }
                    lo = 0, hi = B - 1;
                    while (lo <= hi) {
                        int mid = lo + hi >> 1;
                        ll res = bs[l / B].mx[mid];
                        if (bs[l / B].laz_same != -1) {
                            res = bs[l / B].laz_same;
                        }
                        if (bs[l / B].laz_add != 0) {
                            res += bs[l / B].laz_add;
                        }
                        if (res <= mx) {
                            mx_pos = mid;
                            lo = mid + 1;
                        } else {
                            hi = mid - 1;
                        }
                    }
                    if (mx_pos == mn_pos && mx_pos == -1) {
                        if (bs[l / B].laz_same != -1) {
                            ll val = bs[l / B].laz_same + bs[l / B].laz_add;
                            ans += val * val * B;
                        } else {
                            ans += bs[l / B].pre_val[B - 1];
                            ans += (bs[l / B].pre_mn[B - 1] + bs[l / B].pre_mx[B - 1]) * bs[l / B].laz_add;
                            ans += bs[l / B].laz_add * bs[l / B].laz_add * B;
                        }
                    } else if (mn_pos <= mx_pos) {
                        ans += mn * mx * (mn_pos + 1);
                        ll lhs = bs[l / B].pre_mn[mx_pos];
                        if (mn_pos != -1) lhs -= bs[l / B].pre_mn[mn_pos];
                        if (bs[l / B].laz_same != -1) {
                            lhs = bs[l / B].laz_same * (mx_pos - mn_pos);
                        }
                        if (bs[l / B].laz_add != 0) {
                            lhs += (mx_pos - mn_pos) * bs[l / B].laz_add;
                        }
                        ans += lhs * mx;
                        lhs = bs[l / B].pre_val[B - 1];
                        if (mx_pos != -1) lhs -= bs[l / B].pre_val[mx_pos];
                        if (bs[l / B].laz_same != -1) {
                            ll val = bs[l / B].laz_same;
                            if (bs[l / B].laz_add != 0) val += bs[l / B].laz_add;
                            lhs = val * val * (B - 1 - mx_pos);
                        } else {
                            if (bs[l / B].laz_add != 0) {
                                lhs += bs[l / B].laz_add * bs[l / B].laz_add * (B - 1 - mx_pos);
                                lhs += bs[l / B].laz_add * (bs[l / B].pre_mn[B - 1] + bs[l / B].pre_mx[B - 1]);
                                if (mx_pos != -1) lhs -= bs[l / B].laz_add * (bs[l / B].pre_mn[mx_pos] + bs[l / B].pre_mx[mx_pos]);
                            } else {

                            }
                        }
                        ans += lhs;
                    } else {
                        ans += mn * mx * (mx_pos + 1);
                        ll lhs = bs[l / B].pre_mx[mn_pos];
                        if (mx_pos != -1) lhs -= bs[l / B].pre_mx[mx_pos];
                        if (bs[l / B].laz_same != -1) {
                            lhs = bs[l / B].laz_same * (mn_pos - mx_pos);
                        }
                        if (bs[l / B].laz_add != 0) {
                            lhs += (mn_pos - mx_pos) * bs[l / B].laz_add;
                        }
                        ans += lhs * mn;
                        lhs = bs[l / B].pre_val[B - 1];
                        if (mn_pos != -1) lhs -= bs[l / B].pre_val[mn_pos];
                        if (bs[l / B].laz_same != -1) {
                            ll val = bs[l / B].laz_same;
                            if (bs[l / B].laz_add != 0) val += bs[l / B].laz_add;
                            lhs = val * val * (B - 1 - mn_pos);
                        } else {
                            if (bs[l / B].laz_add != 0) {
                                lhs += bs[l / B].laz_add * bs[l / B].laz_add * (B - 1 - mn_pos);
                                lhs += bs[l / B].laz_add * (bs[l / B].pre_mn[B - 1] + bs[l / B].pre_mx[B - 1]);
                                if (mn_pos != -1) lhs -= bs[l / B].laz_add * (bs[l / B].pre_mn[mn_pos] + bs[l / B].pre_mx[mn_pos]);
                            } else {

                            }
                        }
                        ans += lhs;
                    }
                    ll val = bs[l / B].mx[B - 1];
                    if (bs[l / B].laz_same != -1) {
                        val = bs[l / B].laz_same;
                    }
                    if (bs[l / B].laz_add != 0) {
                        val += bs[l / B].laz_add;
                    }
                    mx = max(mx, val);
                    val = bs[l / B].mn[B - 1];
                    if (bs[l / B].laz_same != -1) {
                        val = bs[l / B].laz_same;
                    }
                    if (bs[l / B].laz_add != 0) {
                        val += bs[l / B].laz_add;
                    }
                    mn = min(mn, val);
                    l += B;
                }
                for (int i = l;i<=r;i++) {
                    mn = min(mn, bs[i / B].ord[i % B]);
                    mx = max(mx, bs[i / B].ord[i % B]);
                    ans += mn * mx;
                }
            }
            cout << ans << '\n';
        }
    }
    return;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    solve();
    return 0;
}
