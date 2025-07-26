```cpp
#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db double
#define ldb long double
#define max(A, B) (A > B ? A : B)
#define min(A, B) (A < B ? A : B)
#define pir(X) pair<X, X>
#define mpr(A, B) make_pair(A, B)
#define fr first
#define sc second
#define sq(x) ((x) * (x))
using namespace std;

template <class T1, class T2>
void cmax(T1 &a, T2 b) { a = max(a, b); }

template <class T1, class T2>
void cmin(T1 &a, T2 b) { a = min(a, b); }

const int N=1e5+10;
int n,m, a[N], b[N], pls[N];
struct segT{
#define root 1,1,n
#define lson (pos<<1)
#define rson (lson|1)
#define mid ((l+r)>>1)
#define l_son lson,l,mid
#define r_son rson,mid+1,r
#define This pos,l,r
	int data[N<<2];
	void upd(int pos){
		data[pos] = max(data[lson], data[rson]);
	}
	void update(int pos,int l,int r,int ip,int val){
		if(l==r){
			assert(ip==l);
			data[pos] += val;
			return;
		}
		if(ip<=mid) update(l_son,ip,val);
		if(ip>=mid+1) update(r_son,ip,val);
		upd(pos);
	}
	int qry(int pos,int l,int r,int L,int R){
		if(L<=l&&r<=R) return data[pos];
		if(R<=mid) return qry(l_son,L,R);
		if(L>=mid+1) return qry(r_son,L,R);
		return max(qry(l_son,L,R), qry(r_son,L,R));
	}
	void build(int pos,int l,int r){
		if(l==r){
			data[pos] = b[l];
			return;
		}
		build(l_son),build(r_son),upd(pos);
	}
}segT;
struct TREE{
	vector<int> t[N];
	int siz[N], son[N], top[N], dfn[N], fa[N], dep[N], dfs_cnt;
	void dfs1(int u,int pre){
		siz[u] = 1, dep[u] = dep[pre]+1;
		fa[u] = pre;
		for(int v:t[u]) if(v!=pre){
			dfs1(v, u);
			siz[u] += siz[v];
			if(siz[v] > siz[son[u]]) son[u] = v;
		}
	}
	void dfs2(int u,int topf){
		top[u] = topf, dfn[u] = ++dfs_cnt;
		if(son[u]==0) return;
		dfs2(son[u], topf);
		
		for(int v:t[u]) 
			if(v!=fa[u] && v!=son[u]) dfs2(v, v);
	}
	void init(){
		for(int i=1;i<=n;i++){
			t[i].clear();
			son[i] = 0;
		}
		for(int i=1;i<n;i++){
			int x, y;
			cin>>x>>y;
			t[x].push_back(y),t[y].push_back(x);
		}
		dfs_cnt = 0, dep[1] = 0;
		dfs1(1, 1);
		dfs2(1, 1);
		for(int i=1;i<=n;i++) b[dfn[i]] = a[i];
		segT.build(root);
	}
	int GET_ANS(int x,int y){
		int ret = 0;
		auto upd=[&](int x) -> void {
			if(x!=fa[x] && pls[fa[x]]>0) cmax(ret, segT.qry(root,dfn[x],dfn[x])+pls[fa[x]]);	
		};
		
		while(top[x] != top[y]){
			if(dep[top[x]] < dep[top[y]]) swap(x, y);
			cmax(ret, segT.qry(root, dfn[top[x]], dfn[x]));
			upd(top[x]);
			x=fa[top[x]];
		}
		if(dep[x] > dep[y]) swap(x,y);
		
		if(x==top[x]) upd(x);
		cmax(ret, segT.qry(root,dfn[x],dfn[y]));
		
	
		return ret;
	}
	void update(int x,int v){
		pls[x] += v;
		if(x!=fa[x]) segT.update(root, dfn[fa[x]], v);
		if(son[x]!=0) segT.update(root, dfn[son[x]], v);
	}
}sptT;

void slv()
{
	cin>>n>>m;
	for(int i=1;i<=n;i++) cin>>a[i], pls[i]=0;
	sptT.init();
	while(m--){
		int opt,x,y;
		cin>>opt>>x>>y;
		if(opt==1){
			cout<<sptT.GET_ANS(x,y)<<'\n';
		}
		if(opt==2){
			sptT.update(x,y);
		}
	}
	return;
}

void fst_IO()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);
}

int main()
{
	fst_IO();
	int T;
	cin >> T;
	while (T--)
	{
		slv();
	}
	return 0;
}
```

copyed from hdu summer3 t11009
