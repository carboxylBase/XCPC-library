// 【框架】Info/Node、合并、标记及业务查询需按题填写；空接口不是可直接使用的实现。
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll N = 2000000;

/*
一个易错点， 树状数组上二分的时候， 注意那个， 这个树状数组的下标是从 1 开始的
这很重要！！！ 离散化的时候要对齐了
*/
struct BIT {
    int n;
    struct Info {
        Info() {
            
        }
        Info operator+(const Info& A) const {
            throw logic_error("TODO: BIT Info merge");
        }
        Info operator-(const Info& A) const {
            throw logic_error("TODO: BIT Info subtract");
        } 
    };
    vector<Info> infos;
    void init(int N_) {
        n = N_ + 1;
        infos.assign(n + 1, Info());
    }
    int lowbit(int x) {return x & -x;}
    void update(int x, ll v) {
        throw logic_error("TODO: BIT update");
        x ++;
        while (x <= n) {

            x += lowbit(x);
        }
    }
    Info query(int x) {
        x ++;
        Info res;
        while (x) {
            res = res + infos[x];
            x -= lowbit(x); 
        }
        return res;
    }
    Info query(int l, int r) {
        if (l > r) return Info();
        return query(r) - query(l - 1);
    }
};
