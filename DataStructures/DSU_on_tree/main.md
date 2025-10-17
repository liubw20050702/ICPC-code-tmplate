## dsu on tree

[彩色的树 ： dsu + stl乱搞](https://codeforces.com/gym/104081/problem/L)

```cpp
const int N = 1e5 + 10;
multiset<int> st;
vector<int> t[N];

int n, m, k, c[N], ans[N];

int dep[N], siz[N], big_son[N];

int cnt = 0;
priority_queue<pir(int)> q;
void ADD(int dep, int col){
	q.push(mpr(dep, col));
	if(st.find(col) == st.end()) cnt++;
	st.insert(col);
}
void DEL(int col){
	st.erase(st.find(col));
	if(st.find(col) == st.end()) cnt--;
}
void reset(int dep){
	while(!q.empty()){
		auto [cur_dep, col] = q.top();
		if(cur_dep - dep <= k) break;
		else{
			q.pop();
			DEL(col);
		}
	}
}
void clear(){
	while(!q.empty()){
		auto [cur_dep, col] = q.top();
		q.pop();
		DEL(col);
	}
}

int dfn[N], undfn[N], dfn_cnt = 0;
void dfs1(int u, int fa){
	siz[u] = 1, dep[u] = dep[fa] + 1;
	dfn[u] = ++dfn_cnt, undfn[dfn_cnt] = u;
	
	for(int v: t[u]) if(v!=fa){
		dfs1(v, u);
		if(big_son[u] == 0 || siz[big_son[u]] < siz[v]) big_son[u] = v;
		siz[u] += siz[v];
	}
}
void dfs2(int u, int fa, bool keep){
	for(int v: t[u]) if(v!=fa && v != big_son[u]) dfs2(v, u, false);
	if(big_son[u]) dfs2(big_son[u], u, true);
	
	for(int v: t[u]) if(v!=fa && v!= big_son[u]){
		for(int i = dfn[v]; i < dfn[v] + siz[v]; i++){
			int pos = undfn[i];
			if(dep[pos] - dep[u] <= k) ADD(dep[pos], c[pos]);
		}
	}
	
	ADD(dep[u], c[u]);
	reset(dep[u]);
	ans[u] = cnt;
	if(keep == false) clear();
}

int main()
{
	fst_IO();
	cin>>n>>k;
	for(int i=1;i<=n;i++) cin>>c[i];
	for(int i=1;i<n;i++){
		int x, y;
		cin>>x>>y;
		t[x].push_back(y), t[y].push_back(x);
	}
	dfs1(1, 1);
	dfs2(1, 1, true);
	
	cin>>m;
	while(m--){
		int x;
		cin>>x;
		cout<<ans[x]<<'\n';
	}
	```
	return 0;
}
```
