整体二分模板 https://www.luogu.com.cn/problem/P3834

离线询问，区间求第 &k& 小

```
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

struct Data{
	int p, x;
};
struct qry{
	int l, r, k, id;
};
const int N = 2e5+10;
int n, m, ans[N];

struct ary_t{
#define lb(x) (x&-x)
	int t[N];
	void upd(int x, int v){
		while(x<=n){
			t[x] += v;
			x += lb(x);
		}
	}
	int qry(int x){
		int ret = 0;
		while(x){
			ret += t[x];
			x -= lb(x);
		}
		return ret;
	}
}T;

void slove(int L, int R, vector<Data> a, vector<qry> q){
	int mid = (L + R)>>1;
	if(L == R){
		for(qry cur : q) ans[cur.id] = L;
		return;
	}
	
	vector<Data> al, ar;
	vector<qry> ql, qr;
	for(auto [p, x] : a){
		if(x <= mid) al.push_back({p, x}), T.upd(p, 1);
		else ar.push_back({p, x});
	}
	
	for(auto [l, r, k, id] : q){
		int num = T.qry(r) - T.qry(l-1);

		if (num < k) qr.push_back({l, r, k-num, id});
		else ql.push_back({l, r, k , id});
	}
	
	for(auto [p, x] : al) T.upd(p, -1);
	slove(L, mid, al, ql), slove(mid+1, R, ar, qr);
	
	return;
}

void slv()
{
	cin>>n>>m;
	vector<Data> a(n);
	vector<int> aa(n);
	for(int i=0;i<n;i++){
		cin>>a[i].x;
		aa[i] = a[i].x;
		a[i].p = i+1;
	}
	int cnt = 0;
	sort(aa.begin(), aa.end());
	map<int,int> mp, de_mp;
	for(int i=0;i<n;i++){
		if(!mp.count(aa[i])) mp[aa[i]] = ++cnt, de_mp[cnt] = aa[i];
	}
	for(int i=0;i<n;i++) a[i].x = mp[a[i].x];
	
	vector<qry> q(m);
	for(int i=0;i<m;i++){
		int l, r, k;
		cin>>l>>r>>k;
		q[i] = {l, r, k, i+1};
	}
	
	slove(1, cnt, a, q);
	for(int i=1;i<=m;i++) cout<<de_mp[ans[i]]<<'\n';
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
	int T = 1;
	while (T--)
	{
		slv();
	}
	return 0;
}
```
