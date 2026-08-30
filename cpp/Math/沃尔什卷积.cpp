#include <vector>
#include <algorithm>

using namespace std;

template <typename T>
struct FWTSolver {
    enum Type { OR = 0, AND = 1, XOR = 2 };

    // 核心变换
    static void transform(vector<T>& a, Type type, bool inv) {
        int n = a.size();
        for (int h = 1; h < n; h <<= 1) {
            for (int j = 0; j < n; j += (h << 1)) {
                for (int k = 0; k < h; k++) {
                    T &x = a[j + k], &y = a[j + k + h];
                    if (type == OR) {
                        if (!inv) y += x; else y -= x;
                    } else if (type == AND) {
                        if (!inv) x += y; else x -= y;
                    } else { // XOR
                        T u = x, v = y;
                        x = u + v; y = u - v;
                        if (inv) { 
                            // 注意：如果是整数或浮点数直接 /2
                            // 如果是 ModInt，此处应乘上 2 的逆元
                            x /= 2; y /= 2; 
                        }
                    }
                }
            }
        }
    }

    // 卷积通用接口
    static vector<T> convolution(vector<T> a, vector<T> b, Type type) {
        int n = 1, sz = max(a.size(), b.size());
        while (n < sz) n <<= 1;
        a.resize(n, T(0)); b.resize(n, T(0));

        transform(a, type, false);
        transform(b, type, false);
        for (int i = 0; i < n; i++) a[i] *= b[i];
        transform(a, type, true);
        return a;
    }

    // 快捷调用
    static vector<T> convOR(vector<T> a, vector<T> b) { return convolution(a, b, OR); }
    static vector<T> convAND(vector<T> a, vector<T> b) { return convolution(a, b, AND); }
    static vector<T> convXOR(vector<T> a, vector<T> b) { return convolution(a, b, XOR); }
};