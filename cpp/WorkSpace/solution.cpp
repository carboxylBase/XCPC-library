#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const ll MOD = 998244353;
ll MMOD(ll x){x%=MOD; if(x<0) x+=MOD; return x;}
ll qpow(ll base,ll k){ ll res=1;base%=MOD; while(k){ if(k&1) res=res*base%MOD; base=base*base%MOD; k>>=1;} return res; }
struct Complex3{ ll x,y; Complex3(ll X_=0,ll Y_=0){x=MMOD(X_); y=MMOD(Y_);} Complex3 operator+(const Complex3& o)const{return Complex3(MMOD(x+o.x),MMOD(y+o.y));} Complex3 operator-(const Complex3& o)const{return Complex3(MMOD(x-o.x),MMOD(y-o.y));} Complex3 operator*(const Complex3& o)const{ ll a=x,b=y,c=o.x,d=o.y; ll rx = MMOD((a*c%MOD - b*d%MOD)); ll ry = MMOD((a*d%MOD + b*c%MOD - b*d%MOD)); return Complex3(rx,ry);} bool operator==(const Complex3& o) const { return x==o.x && y==o.y; } Complex3 inv() const { ll a=x,b=y; ll n = MMOD((a*a%MOD - a*b%MOD + b*b%MOD)); n = qpow(n, MOD-2); return Complex3(MMOD((a - b)%MOD * n %MOD), MMOD((-b)%MOD * n %MOD)); } Complex3 operator/(const Complex3& o) const { return (*this) * o.inv(); } };
Complex3 zeroC(){ return Complex3(0,0); }

Complex3 det_gauss(vector<vector<Complex3>> A){ int n = (int)A.size(); int sign = 0; for(int i=0;i<n;i++){ int piv = i; for(int j=i;j<n;j++) if(!(A[j][i]==zeroC())){ piv=j; break; } if(A[piv][i]==zeroC()){ return zeroC(); } if(piv!=i){ swap(A[piv], A[i]); sign ^= 1; } Complex3 inv = A[i][i].inv(); for(int j=i+1;j<n;j++){ if(A[j][i]==zeroC()) continue; Complex3 factor = A[j][i] / A[i][i]; for(int k=i;k<n;k++) A[j][k] = A[j][k] - factor * A[i][k]; } } Complex3 res = Complex3(1,0); for(int i=0;i<n;i++) res = res * A[i][i]; if(sign) res = Complex3(MMOD(-res.x), MMOD(-res.y)); return res; }

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int n,m; if(!(cin>>n>>m)) return 0;
    struct E{int u,v; ll w; int t;};
    vector<E> es(m);
    for(int i=0;i<m;i++) cin>>es[i].u>>es[i].v>>es[i].w>>es[i].t;
    auto build_laplacian = [&](int power)->vector<vector<Complex3>>{
        vector<vector<Complex3>> L(n, vector<Complex3>(n, zeroC()));
        for(auto &e: es){
            int exp = (e.t==0?1:(e.t==1?-1:0));
            int k = ((power * exp) % 3 + 3) % 3;
            Complex3 val;
            if(k==0) val = Complex3(e.w,0);
            else if(k==1) val = Complex3(0,e.w);
            else val = Complex3(MMOD(-e.w), MMOD(-e.w));
            int u = e.u-1, v = e.v-1;
            L[u][u] = L[u][u] + val;
            L[v][v] = L[v][v] + val;
            L[u][v] = L[u][v] - val;
            L[v][u] = L[v][u] - val;
        }
        // return minor (n-1)x(n-1)
        vector<vector<Complex3>> M(n-1, vector<Complex3>(n-1, zeroC()));
        for(int i=0;i<n-1;i++) for(int j=0;j<n-1;j++) M[i][j]=L[i][j];
        return M;
    };
    Complex3 sum = zeroC();
    for(int p=0;p<3;p++){
        auto M = build_laplacian(p);
        sum = sum + det_gauss(M);
    }
    ll inv3 = qpow(3, MOD-2);
    // divide both components by 3 (imag should be 0)
    sum.x = sum.x * inv3 % MOD;
    sum.y = sum.y * inv3 % MOD;
    cout << sum.x << "\n";
    return 0;
}
