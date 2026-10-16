from pathlib import Path
import subprocess,concurrent.futures,json
import tempfile, shutil, os
root=Path(__file__).resolve().parents[1]/'cpp'
temp=tempfile.TemporaryDirectory(prefix='xcpc-template-tests-')
base=Path(temp.name)
cc=os.environ.get('CXX') or shutil.which('g++')
if not cc: raise SystemExit('Set CXX to a GCC-compatible C++ compiler path.')
def read(f):
 return (root/f).read_text(encoding='utf-8-sig').replace('#include "all.hpp"','#include <bits/stdc++.h>')
def scan(p):
 f=str(p.relative_to(root));dest=base/('check_'+p.parent.name+'_'+p.name);dest.write_text(read(f),encoding='utf-8')
 r=subprocess.run([cc,'-std=c++17','-fsyntax-only',str(dest)],capture_output=True,encoding='utf-8',errors='replace')
 return f,r.returncode,r.stderr
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:r=list(pool.map(scan,[p for p in root.rglob('*.cpp') if 'WorkSpace' not in p.parts]))
print('syntax:',len(r),'failed:',sum(x[1]!=0 for x in r),flush=True)
for f,c,e in r:
 if c:print(f,e[:1200],flush=True)
assert all(c==0 for _,c,_ in r)
tests={
'lca':('Graph/Lca.cpp','vector<vector<int>> g(3);g[1]={2};g[2]={1};solverLcaBinary.init(g,1);assert(solverLcaBinary.lca(1,2)==1);'),
'trie':('String/字典树.cpp','Trie t;mt19937 rng(42);vector<ll> a;for(int i=0;i<500;i++){ll v=rng()%100000;a.push_back(v);t.insert(v);ll x=rng()%100000,z=0;for(ll y:a)z=max(z,x^y);assert(t.getMx(x)==z);}'),
'basis':('Math/线性基.cpp','LinearBasis b;assert(b.get_mn()==-1);b.insert(1);assert(!b.flag0);b.insert(1);assert(b.flag0&&b.get_mn()==0);b.init();b.insert(0);assert(b.flag0);'),
'floor':('Math/类欧几里得算法.cpp','for(int n=0;n<10;n++)for(int m=1;m<10;m++)for(int a=-20;a<=20;a++)for(int b=-20;b<=20;b++){long long z=0;for(int i=0;i<n;i++){int x=a*i+b;z+=x>=0?x/m:-((-x+m-1)/m);}assert(floor_sum_mod(n,m,a,b,97)==(z%97+97)%97);}'),
'lichao':('DataStructure/李超树.cpp','static LiChaoTree t;t.init(100);mt19937 rng(4);vector<array<int,3>> a;for(int i=0;i<200;i++){int l=rng()%100+1,r=rng()%100+1;if(l>r)swap(l,r);int y=int(rng()%200)-100;a.push_back({l,r,y});t.add(l,y,r,y);t.update(1,1,100,l,r,i+1);for(int x=1;x<=100;x++){double z=-1e18;for(auto v:a)if(v[0]<=x&&x<=v[1])z=max(z,double(v[2]));assert(t.query(1,1,100,x).first==z);}}'),
'flow':('Graph/流.cpp','MinCostMaxFlow<ll> f(2,10);f.addEdge(1,2,1,2000000000LL);assert(f.run(1,2)==make_pair(1LL,2000000000LL));f.init(4);f.addEdge(1,4,2,3);assert(f.run(1,4)==make_pair(2LL,6LL));'),
'pam':('String/回文自动机.cpp','mt19937 rng(2);for(int z=0;z<200;z++){PAM::s="";for(int i=0;i<20;i++)PAM::s+=char(\'a\'+rng()%3);PAM::init();set<string>a;for(int i=0;i<20;i++)for(int j=i;j<20;j++){auto s=PAM::s.substr(i,j-i+1),t=s;reverse(t.begin(),t.end());if(s==t)a.insert(s);}assert(PAM::tot-1==a.size());}'),
'sam':('String/后缀自动机(SAM).cpp','mt19937 rng(3);for(int z=0;z<200;z++){string s;for(int i=0;i<20;i++)s+=char(\'a\'+rng()%3);SAM::sam_init();SAMLowercase::sam_init();for(char c:s){SAM::sam_extend(c);SAMLowercase::sam_extend(c);}set<string>a;for(int i=0;i<20;i++)for(int j=i;j<20;j++)a.insert(s.substr(i,j-i+1));long long x=0,y=0;for(int i=1;i<SAM::sz;i++)x+=SAM::st[i].len-SAM::st[SAM::st[i].link].len;for(int i=1;i<SAMLowercase::sz;i++)y+=SAMLowercase::st[i].len-SAMLowercase::st[SAMLowercase::st[i].link].len;assert(x==a.size()&&y==a.size());}'),
'splay':('DataStructure/伸展树.cpp','static Splay<int> t;t.init();t.insert(5);t.insert(5);assert(t.nodes[t.rt].siz==4);'),
'dijkstra':('Graph/最短路.cpp','static Graph<ll> g;static Dijkstra<ll>d;g.init(5,5);g.add(1,2,1);g.add(2,3,2147483647LL);g.add(1,4,2);g.add(4,3,1);g.add(3,5,1);d.g=&g;d.dijkstra(1);assert(d.dist[5]==4);'),
'matrix':('Math/矩阵.cpp','Mat<ll>a(2,2);a.a={{1,1},{1,0}};a=a;assert(a.a[0][1]==1);assert(qpow(a,0).a==Mat<ll>::identity(2).a);auto b=Mat<ll>::identity(2);for(int i=0;i<100;i++){assert(qpow(a,i).a==b.a);b*=a;}'),
'closest':('Geometry/平面最近点对.cpp','vector<pair<int,int>> a={{-1500000000,0},{1500000000,0}};assert(closest_pair_sq(a)==9000000000000000000LL);a.assign(200000,{1,1});assert(closest_pair_sq(a)==0);mt19937 rng(7);for(int z=0;z<200;z++){a.clear();for(int i=0;i<50;i++)a.push_back({int(rng()%1000),int(rng()%1000)});long long best=LLONG_MAX;for(int i=0;i<50;i++)for(int j=0;j<i;j++){long long x=a[i].first-a[j].first,y=a[i].second-a[j].second;best=min(best,x*x+y*y);}assert(closest_pair_sq(a)==best);}'),
'hld':('Graph/树链剖分.cpp','static TreeChainSeg t;vector<vector<int>>g(2);t.init(g);assert(t.dfn[1]==1);'),
'heap':('DataStructure/左偏树.cpp','static LeftLeaningHeap h;int rt=0;for(int i=1;i<=1000;i++){h.nodes[i].fa=i;h.nodes[i].val=i;rt=h.merge(rt,i);}auto dfs=[&](auto&&self,int u)->int{if(!u)return -1;int l=self(self,h.nodes[u].ls),r=self(self,h.nodes[u].rs);assert(l>=r&&h.nodes[u].dist==r+1);return r+1;};dfs(dfs,rt);'),
'cartesian':('DataStructure/笛卡尔树.cpp','CartesianTree<int> t;vector<int>a={0,-2,-3,1};t.buildMinHeap(a);assert(t.nodes[0].ch[1]==2&&t.nodes[2].ch[0]==1&&t.nodes[2].ch[1]==3);'),
}
for name,(f,body) in tests.items():
 src=read(f)
 if name=='splay':src=src.split('struct Node {\n    ll val;')[0]
 p=base/('fixedtest_'+name+'.cpp');p.write_text(src+'\nint main(){'+body+'}\n',encoding='utf-8');exe=p.with_suffix('.exe')
 r=subprocess.run([cc,'-std=c++17','-O2',str(p),'-o',str(exe)],capture_output=True,encoding='utf-8',errors='replace');assert r.returncode==0,(name,r.stderr)
 r=subprocess.run([str(exe)],capture_output=True,timeout=10);assert r.returncode==0,(name,r.stderr)
 print('PASS',name,flush=True)
