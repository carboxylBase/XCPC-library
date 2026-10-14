#include <bits/stdc++.h>
using namespace std;

// 三种 build_sa 方法任选一种，然后调用 build_height。
// s 为 0-based；sa/rk/height 为 1-based，所有方法均可重复调用。
// build_sa()        : 倍增 + sort，O(n log^2 n)
// build_sa_nlogn()  : 倍增 + 计数排序，O(n log n)
// build_sa_linear() : SA-IS，字节字符集下 O(n)
// build_height()    : Kasai，O(n)
// SA-IS 核心改编自 AtCoder Library (CC0):
// https://github.com/atcoder/ac-library/blob/master/atcoder/string.hpp
struct SA {
    string s;
    int n = 0;
    vector<int> sa, rk, height;
    void build_sa() {
        n = s.size();
        sa.assign(n + 1, 0);
        rk.assign(2 * n + 2, 0);
        height.assign(n + 1, 0);
        vector<int> oldrk(2 * n + 2, 0);

        for (int i = 1; i <= n; i++) {
            sa[i] = i;
            rk[i] = (unsigned char)s[i - 1] + 1;
        }

        for (int w = 1; w < n; w <<= 1) {
            sort(sa.begin() + 1, sa.end(), [&](int x, int y) {
                if (rk[x] != rk[y]) return rk[x] < rk[y];
                return rk[x + w] < rk[y + w];
            });

            oldrk = rk;
            int p = 0;
            for (int i = 1; i <= n; i++) {
                int x = sa[i], y = sa[i - 1];
                if (i == 1 || oldrk[x] != oldrk[y] ||
                    oldrk[x + w] != oldrk[y + w]) {
                    ++p;
                }
                rk[x] = p;
            }
            if (p == n) break;
        }

        // 保证 rk 为排名，也处理 n = 1 的情况
        for (int i = 1; i <= n; i++) rk[sa[i]] = i;
    }

    void build_sa_nlogn() {
        n = s.size();
        sa.assign(n + 1, 0);
        rk.assign(2 * n + 2, 0);
        height.assign(n + 1, 0);
        if (!n) return;

        vector<int> oldrk(2 * n + 2, 0), id(n + 1);
        vector<int> cnt(max(n, 256) + 1, 0);
        int m = 256;
        for (int i = 1; i <= n; i++) {
            rk[i] = (unsigned char)s[i - 1] + 1;
            ++cnt[rk[i]];
        }
        for (int i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
        for (int i = n; i >= 1; i--) sa[cnt[rk[i]]--] = i;

        for (int w = 1; w < n; w <<= 1) {
            // id 按第二关键字 rk[i+w] 排序：
            // 越界的第二段排名为 0；其余由上一轮 sa 整体左移 w 得到。
            int p = 0;
            for (int i = n - w + 1; i <= n; i++) id[++p] = i;
            for (int i = 1; i <= n; i++) {
                if (sa[i] > w) id[++p] = sa[i] - w;
            }

            // 按第一关键字稳定计数排序，保留第二关键字顺序。
            fill(cnt.begin(), cnt.begin() + m + 1, 0);
            for (int i = 1; i <= n; i++) ++cnt[rk[i]];
            for (int i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
            for (int i = n; i >= 1; i--) {
                sa[cnt[rk[id[i]]]--] = id[i];
            }

            oldrk = rk;
            p = 0;
            for (int i = 1; i <= n; i++) {
                int x = sa[i], y = sa[i - 1];
                if (i == 1 || oldrk[x] != oldrk[y] ||
                    oldrk[x + w] != oldrk[y + w]) {
                    ++p;
                }
                rk[x] = p;
            }
            m = p;
            if (m == n) break;
        }
        for (int i = 1; i <= n; i++) rk[sa[i]] = i;
    }

    void build_sa_linear() {
        n = s.size();
        vector<int> a(n);
        for (int i = 0; i < n; i++) a[i] = (unsigned char)s[i];
        vector<int> order = sa_is(a, 255);
        sa.assign(n + 1, 0);
        rk.assign(n + 1, 0);
        height.assign(n + 1, 0);
        for (int i = 1; i <= n; i++) {
            sa[i] = order[i - 1] + 1;
            rk[sa[i]] = i;
        }
    }

    // height[i] = LCP(sa[i - 1], sa[i])，height[1] = 0
    void build_height() {
        height.assign(n + 1, 0);
        for (int i = 1, h = 0; i <= n; i++) {
            if (rk[i] == 1) {
                h = 0;
                continue;
            }
            if (h) --h;
            int j = sa[rk[i] - 1];
            while (i + h <= n && j + h <= n &&
                   s[i + h - 1] == s[j + h - 1]) {
                ++h;
            }
            height[rk[i]] = h;
        }
    }

private:
    // ls[i] 为 S 型标记；LMS 是前一位 L 型、当前位置 S 型的位置。
    // 先按 LMS 排序，再诱导 L 型和 S 型后缀；LMS 子串命名后递归。
    vector<int> sa_is(const vector<int>& s, int upper) {
        int n = int(s.size());
        if (n == 0) return {};
        if (n == 1) return {0};
        if (n == 2) {
            if (s[0] < s[1]) {
                return {0, 1};
            } else {
                return {1, 0};
            }
        }
        vector<int> sa(n);
        vector<bool> ls(n);
        for (int i = n - 2; i >= 0; i--) {
            ls[i] = (s[i] == s[i + 1]) ? ls[i + 1] : (s[i] < s[i + 1]);
        }
        vector<int> sum_l(upper + 1), sum_s(upper + 1);
        for (int i = 0; i < n; i++) {
            if (!ls[i]) {
                sum_s[s[i]]++;
            } else {
                sum_l[s[i] + 1]++;
            }
        }
        for (int i = 0; i <= upper; i++) {
            sum_s[i] += sum_l[i];
            if (i < upper) sum_l[i + 1] += sum_s[i];
        }
    
        auto induce = [&](const vector<int>& lms) {
            fill(sa.begin(), sa.end(), -1);
            vector<int> buf(upper + 1);
            copy(sum_s.begin(), sum_s.end(), buf.begin());
            for (auto d : lms) {
                if (d == n) continue;
                sa[buf[s[d]]++] = d;
            }
            copy(sum_l.begin(), sum_l.end(), buf.begin());
            sa[buf[s[n - 1]]++] = n - 1;
            for (int i = 0; i < n; i++) {
                int v = sa[i];
                if (v >= 1 && !ls[v - 1]) {
                    sa[buf[s[v - 1]]++] = v - 1;
                }
            }
            copy(sum_l.begin(), sum_l.end(), buf.begin());
            for (int i = n - 1; i >= 0; i--) {
                int v = sa[i];
                if (v >= 1 && ls[v - 1]) {
                    sa[--buf[s[v - 1] + 1]] = v - 1;
                }
            }
        };
    
        vector<int> lms_map(n + 1, -1);
        int m = 0;
        for (int i = 1; i < n; i++) {
            if (!ls[i - 1] && ls[i]) {
                lms_map[i] = m++;
            }
        }
        vector<int> lms;
        lms.reserve(m);
        for (int i = 1; i < n; i++) {
            if (!ls[i - 1] && ls[i]) {
                lms.push_back(i);
            }
        }
    
        induce(lms);
    
        if (m) {
            vector<int> sorted_lms;
            sorted_lms.reserve(m);
            for (int v : sa) {
                if (lms_map[v] != -1) sorted_lms.push_back(v);
            }
            vector<int> rec_s(m);
            int rec_upper = 0;
            rec_s[lms_map[sorted_lms[0]]] = 0;
            for (int i = 1; i < m; i++) {
                int l = sorted_lms[i - 1], r = sorted_lms[i];
                int end_l = (lms_map[l] + 1 < m) ? lms[lms_map[l] + 1] : n;
                int end_r = (lms_map[r] + 1 < m) ? lms[lms_map[r] + 1] : n;
                bool same = true;
                if (end_l - l != end_r - r) {
                    same = false;
                } else {
                    while (l < end_l) {
                        if (s[l] != s[r]) {
                            break;
                        }
                        l++;
                        r++;
                    }
                    if (l == n || s[l] != s[r]) same = false;
                }
                if (!same) rec_upper++;
                rec_s[lms_map[sorted_lms[i]]] = rec_upper;
            }
    
            auto rec_sa =
                sa_is(rec_s, rec_upper);
    
            for (int i = 0; i < m; i++) {
                sorted_lms[i] = lms[rec_sa[i]];
            }
            induce(sorted_lms);
        }
        return sa;
    }
} solver;

