#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const ll N = 2000000;
const ll INF = 1e18;

/*
切记 pd 清空旧标记
*/

#define LS (rt << 1)
#define RS (rt << 1 | 1)

struct SegTree{
    struct Node{
        Node() {

        }
    }nodes[N*4];

    Node merge(Node L,Node R){
        Node M;
        return M;
    }
    void build(int rt,int l,int r){
        if (l == r){
            return;
        }
        int mid = l + r >> 1;
        build(LS,l,mid),build(RS,mid+1,r);
        nodes[rt] = merge(nodes[LS],nodes[RS]);
    }
    void pd(int rt){

    }
    void update(int rt,int l,int r,int ql,int qr,int val){
        if (ql > qr || l > qr || ql > r) return;
        if (ql <= l && r <= qr){
            return;
        }
        int mid = l+r>>1;
        pd(rt);
        if (ql <= mid){
            update(LS,l,mid,ql,qr,val);
        }
        if (qr >= mid + 1){
            update(RS,mid+1,r,ql,qr,val);
        }
        nodes[rt] = merge(nodes[LS],nodes[RS]);
        return;
    }
    void modify(int rt,int l,int r,int q,int val){
        if (l > q || r < q) return;
        if (l == r) {
            return;
        }
        int mid = l+r>>1;
        pd(rt);
        if (q <= mid){
            modify(LS,l,mid,q,val);
        }
        if (q >= mid + 1){
            modify(RS,mid+1,r,q,val);
        }
        nodes[rt] = merge(nodes[LS],nodes[RS]);
        return;
    }
    Node query(int rt,int l,int r,int ql,int qr){
        if (ql > qr) {
            return Node();
        }
        if (ql <= l && r <= qr){
            return nodes[rt];
        }
        int mid = l+r>>1;
        pd(rt);
        if (ql > mid){
            return query(RS,mid+1,r,ql,qr);
        }else if (qr < mid + 1){
            return query(LS,l,mid,ql,qr);
        }else{
            return merge(query(LS,l,mid,ql,qr),query(RS,mid+1,r,ql,qr));
        }
    }
};

// 调用前先 init
#define LS nodes[rt].ls
#define RS nodes[rt].rs
struct DynamicSegTree {
    struct Node {
        int ls,rs;
        Node() {
            ls = rs = 0;
        }
    } nodes[N];
    int cntNode;
    void init() {
        cntNode = 1;
        nodes[0] = Node();
    }
    int newNode() {
        nodes[cntNode] = Node();
        return cntNode++;
    }
    Node merge(Node L,Node R) {
        Node M;
        // 补充合并方法
        return M;
    }
    void pushUp(int rt) {
        // 更新父节点
        return;
    }
    void pushDown(int rt) {
        // 下传标记
    }
    void update(int& rt,int l,int r,int ql,int qr,int val) {
        if (ql > qr || l > qr || r < ql) return;
        if (!rt) rt = newNode();
        if (ql <= l && r <= qr) {
            // 更新方法
            return;
        }
        int mid = l + r >> 1;
        update(nodes[rt].ls,l,mid,ql,qr,val);
        update(nodes[rt].rs,mid+1,r,ql,qr,val);
        pushUp(rt);
        return;
    }
    Node query(int& rt,int l,int r,int ql,int qr) {
        if (ql > qr || l > qr || r < ql || !rt) return Node();
        if (ql <= l && r <= qr) {
            return nodes[rt];
        }
        int mid = l + r >> 1;
        pushDown(rt);
        return merge(query(LS,l,mid,ql,qr),query(RS,mid+1,r,ql,qr));
    }
    int treeMerge(int& ls,int& rs,int l,int r) {
        if (!ls || !rs) return ls | rs;
        if (l == r) {
            // 叶节点合并方法
            
            return ls;
        }
        int mid = l + r >> 1;
        nodes[ls].ls = treeMerge(nodes[ls].ls,nodes[rs].ls,l,mid);
        nodes[ls].rs = treeMerge(nodes[ls].rs,nodes[rs].rs,mid+1,r);
        pushUp(ls);
        return ls;
    }
};


struct DynamicSegTree {
	#define LS info[rt].ls
	#define RS info[rt].rs
	struct Info {
        int ls, rs;
        Info() {
            ls = rs = 0;
        }
        Info operator+(const Info& A) const {

        }
    } info[N];
	int cntNode = 0;
	
	void init() {
		cntNode = 0;
	};

	int newNode() {
		++cntNode;
		info[cntNode] = Info();
		return cntNode;
	};

	void pushUp(int rt) {
		return;
	}

	void build(int& rt,int l,int r) {
		rt = newNode();
		if (l == r) {
			return;
		}
		int mid = l + r >> 1;
		build(LS,l,mid); build(RS,mid+1,r);
		pushUp(rt);
		return;
	}

    // 主席树版本
	// void update(int& rt,int l,int r,int ql,int qr) {
	// 	if (ql > qr || qr < l || r < ql) return;
	// 	if (ql <= l && r <= qr) {
	// 		return;
	// 	}
	// 	int mid = l + r >> 1;
	// 	if (ql <= mid) {
	// 		++cntNode;
	// 		info[cntNode] = info[LS];
	// 		LS = cntNode;
	// 		update(LS,l,mid,ql,qr);
	// 	}
	// 	if (qr >= mid + 1) {
	// 		++cntNode;
	// 		info[cntNode] = info[RS];
	// 		RS = cntNode;
	// 		update(RS,mid+1,r,ql,qr);
	// 	}
	// 	pushUp(rt);
	// 	return;
	// }

    void update(int& rt,int l,int r,int ql,int qr) {
		if (ql > qr || qr < l || r < ql) return;
        if (!rt) rt = newNode();
		if (ql <= l && r <= qr) {
			return;
		}
		int mid = l + r >> 1;
		if (ql <= mid) {
			update(LS,l,mid,ql,qr);
		}
		if (qr >= mid + 1) {
			update(RS,mid+1,r,ql,qr);
		}
		pushUp(rt);
		return;
	}

	Info query(int rt,int l,int r,int ql,int qr) {
		if (!rt || ql > qr || r < ql || qr < l) return Info();
		if (ql <= l && r <= qr) {
			return info[rt];
		}
		int mid = l + r >> 1;
		if (qr <= mid) {
			return query(LS,l,mid,ql,qr);
		} else if (ql >= mid + 1) {
			return query(RS,mid+1,r,ql,qr);
		} else {
			return query(LS,l,mid,ql,qr) + query(RS,mid+1,r,ql,qr);
		}
	}
} solver;


/*
线段树合并以及线段树分裂的整合代码
细节还是比较多的, 比如说, 注意这个 split 之后的 pushUp
*/
namespace SegTree {
    struct Node {
        int ch[2];
        ll cnt, sum;
        Node() {
            ch[0] = ch[1] = cnt = sum = 0;
        }
    } nodes[N * 4];
    int tp;

    void init() { tp = 0; }

    int New() {
        ++ tp;
        nodes[tp] = Node();
        return tp;
    }

    #define ls(x) nodes[x].ch[0]
    #define rs(x) nodes[x].ch[1]

    void push_up(int rt) {
        nodes[rt].sum = nodes[rt].cnt = 0;
        if (ls(rt)) {
            nodes[rt].sum += nodes[ls(rt)].sum;
            nodes[rt].cnt += nodes[ls(rt)].cnt;
        }
        if (rs(rt)) {
            nodes[rt].sum += nodes[rs(rt)].sum;
            nodes[rt].cnt += nodes[rs(rt)].cnt;
        }
        return;
    }

    void update(int &rt, int l, int r, int q, ll c) {
        if (!rt) rt = New();
        if (l == r) {
            nodes[rt].sum += c;
            nodes[rt].cnt += c;
            return;
        }
        int mid = l + r >> 1;
        if (q <= mid) update(ls(rt), l, mid, q, c);
        else update(rs(rt), mid + 1, r, q, c);
        push_up(rt);
        return;
    }

    ll query(int rt, int l, int r, int ql, int qr) {
        if (!rt) return 0;
        if (qr < l || ql > r || ql > qr) return 0;
        if (ql <= l && r <= qr) return nodes[rt].sum;
        int mid = l + r >> 1;
        return query(ls(rt), l, mid, ql, qr) + query(rs(rt), mid + 1, r, ql, qr);
    }

    int findK(int rt, int l, int r, ll k) {
        if (!rt || nodes[rt].sum < k) return -1;
        if (l == r) return l;
        int mid = l + r >> 1;
        if (nodes[ls(rt)].sum >= k) {
            return findK(ls(rt), l, mid, k);
        } else {
            k -= nodes[ls(rt)].sum;
            return findK(rs(rt), mid + 1, r, k);
        }
    }

    int merge(int x, int y, int l, int r) {
        if (!x || !y) return x + y;
        if (l == r) {
            nodes[x].cnt += nodes[y].cnt;
            nodes[x].sum += nodes[y].sum;
            return x;
        }
        int mid = l + r >> 1;
        ls(x) = merge(ls(x), ls(y), l, mid);
        rs(x) = merge(rs(x), rs(y), mid + 1, r);
        push_up(x);
        return x;
    }

    void split(int x, int& y, int l, int r, ll k) {
        if (!k) return;
        if (!y) y = New();
        if (l == r) {
            // assert(nodes[x].sum >= k);
            nodes[x].cnt -= k;
            nodes[x].sum -= k;
            nodes[y].cnt += k;
            nodes[y].sum += k;
            return;
        }
        int mid = l + r >> 1;
        if (nodes[ls(x)].sum <= k) {
            k -= nodes[ls(x)].sum;
            nodes[y].ch[0] = ls(x);
            ls(x) = 0;
            split(rs(x), rs(y), mid + 1, r, k);
        } else {
            split(ls(x), ls(y), l, mid, k);
        }
        push_up(x); push_up(y);
        return;
    }
};
using namespace SegTree;