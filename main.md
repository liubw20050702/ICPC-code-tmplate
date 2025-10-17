 - 快读
```cpp
template <class T>void read(T &x){
    x=0;
    char c=getchar(),d='0';
    while(c<'0'||c>'9') d=c,c=getchar();
    while(c>='0'&&c<='9') x=(x<<3)+(x<<1)+c-'0',c=getchar();
    if(d=='-') x=-x;
}
template <class T>void wt(T x){
    if(x/10) wt(x/10);
    putchar(x%10+'0');
}
template <class T>void enter(T x){
    if(x<0) x=-x,putchar('-');
    wt(x),putchar('\n');
}
template <class T>void space(T x){
    if(x<0) x=-x,putchar('-');
    wt(x),putchar(' ');
}

```

 - 主席树
```cpp
const int N = 2e5 + 10;
const int inf = 1e9 + 10;

class persistent_SEGtree {
public:
#define mid ((l + r) >> 1)

    int rt[N];
    persistent_SEGtree() { rt[0] = 0, tot = 0; }
    void update(int &pos, int l, int r, int ip, int v) {
        pos = cpy(pos);
        data[pos] += v;
        if (l == r)
            return;
        if (l <= ip && ip <= mid)
            update(ls[pos], l, mid, ip, v);
        if (mid + 1 <= ip && ip <= r)
            update(rs[pos], mid + 1, r, ip, v);
    }
    int qry(int pos1, int pos2, int l, int r, int k) {
        if (l == r)
            return l;
        int rec = data[ls[pos2]] - data[ls[pos1]];
        if (k <= rec)
            return qry(ls[pos1], ls[pos2], l, mid, k);
        else
            return qry(rs[pos1], rs[pos2], mid + 1, r, k - rec);
    }

private:
    static const int M = (N << 5);
    int tot, data[M], ls[M], rs[M];
    int cpy(int pos) {
        int y = ++tot;
        data[y] = data[pos];
        ls[y] = ls[pos], rs[y] = rs[pos];
        return y;
    }
#undef mid
} T;

void fst_IO() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);
}

int main() {
    fst_IO();
    int n, m;
    cin >> n >> m;
    for (int i = 1; i <= n; i++) {
        int x;
        cin >> x;
        T.rt[i] = T.rt[i - 1];
        T.update(T.rt[i], -inf, inf, x, 1);
    }
    while (m--) {
        int l, r, k;
        cin >> l >> r >> k;
        cout << T.qry(T.rt[l - 1], T.rt[r], -inf, inf, k) << '\n';
    }
    return 0;
}
```

 - 费用流
```cpp
class Dinic {
public:
    void add(int x, int y, int z, int zz, int co) {
        to[++tot] = y, cap[tot] = z, cost[tot] = co, nxt[tot] = head[x], head[x] = tot;
        to[++tot] = x, cap[tot] = zz, cost[tot] = -co, nxt[tot] = head[y], head[y] = tot;
    }
    pair<int,int> dinic() {
        co = 0;
        int ans = 0;
        while (bfs()) ans += dfs(s, INF);
        return make_pair(co,ans);
    }
    Dinic(int S,int T){
        s=S,t=T,tot=1;
        memset(head,0,sizeof head);
        memset(nxt,0,sizeof nxt);
        memset(vis,0,sizeof vis);
        while(!q.empty()) q.pop();
    }
private:
    static const int INF = 0x3f3f3f3f;
    static const int  N = 505;
    static const int M = 100010;

    int head[N], to[M<<1], cap[M<<1], nxt[M<<1], tot = 1, cur[M<<1], cost[M<<1];
    int s, t ,co;
    int dis[N];
    bool vis[N];
    queue<int> q;
    bool bfs() {
        for(int i=0;i<N;i++) dis[i]=INF;
        memcpy(cur, head, sizeof(head));
        q.push(s), dis[s] = 0, vis[s] = 1;
        while (!q.empty()) {
            int u = q.front();
            q.pop(), vis[u] = 0;
            for (int i = head[u]; i; i = nxt[i]) {
                int v = to[i];
                if (cap[i] && dis[v] > dis[u] + cost[i]) {
                    dis[v] = dis[u] + cost[i];
                    if (!vis[v])
                        q.push(v), vis[v] = 1;
                }
            }
        }
        return dis[t] != INF;
    }
    int dfs(int u, int flow) {
        if (u == t)
            return flow;
        vis[u] = 1;
        int ans = 0;
        for (int i = cur[u]; i; i = nxt[i]) {
            cur[u] = i;
            int v = to[i];
            if (!vis[v] && cap[i] && dis[v] == dis[u] + cost[i]) {
                int x = dfs(v, min(cap[i], flow - ans));
                if (x)
                    co += x * cost[i], cap[i] -= x, cap[i ^ 1] += x, ans += x;
            }
            if (ans == flow)
                break;
        }
        if (!ans)
            dis[u] = -1;
        vis[u] = 0;
        return ans;
    }
};


int main() {
    int n,m;
    cin>>n>>m;
    Dinic G(1,n);
    for (int i = 1; i <= m; i++) {
        int a, b, c, d;
        cin >> a >> b >> c >> d;
        G.add(a, b, c, 0, d);
    }
    pair<int,int> p = G.dinic();
    cout << p.second << ' ' << p.first << '\n';
    return 0;
}
```
 - 树剖 + 线段树
```cpp
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
```

 - 莫队
```cpp
int b[N],n,m,k,siz;
ll res[N];
struct quest{
    int l;
    int r;
    int id,pos;
    bool operator <(const quest &b)const{
        if(pos==b.pos) return r<b.r;
        else return pos<b.pos;
    }
}a[N];
void slv(){
    cin>>n>>m;
    siz=sqrt(n);
    for(int i=1;i<=n;i++) cin>>b[i];
    s.clear();
    for(int i=1;i<=m;i++){
        cin>>a[i].l>>a[i].r;
        a[i].id=i;
        a[i].pos=(a[i].l-1)/siz+1;
    }
    sort(a+1,a+1+m);
    int l=1,r=0;
    ll ans=0;
    for(int i=1;i<=m;i++){
        while(l>a[i].l) ins(b[--l]);
        while(r<a[i].r) ins(b[++r]);
        while(l<a[i].l) del(b[l++]);
        while(r>a[i].r) del(b[r--]);
        res[a[i].id]=chk();
    }
    for(int i=1;i<=m;i++) 
        if(res[i]) cout<<"YES\n";
        else cout<<"NO\n";
    return;
}
```
 - 倍增求lca
```cpp
const int N=5e5+10;
const int L=21;
int n,m,s,f[N][L],dep[N];
vector<int> t[N];

void dfs(int pos,int fa){
    f[pos][0]=fa;
    for(int i=0;i+1<L;i++) f[pos][i+1]=f[f[pos][i]][i];
    for(int v:t[pos]) if(v!=fa){
        dep[v]=dep[pos]+1;
        dfs(v,pos);
    }
}
int lca(int a,int b){
    if(dep[a]<dep[b]) swap(a,b);
    int l=dep[a]-dep[b];
    for(int i=L-1;i>=0;i--) if(l&(1<<i)) a=f[a][i];
    if(a==b) return a;
    for(int i=L-1;i>=0;i--){
        if(f[a][i]==f[b][i]) continue;
        a=f[a][i],b=f[b][i];
    }
    return f[a][0];
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(0);cout.tie(0);
    cin>>n>>m>>s;
    for(int i=1;i<n;i++){
        int x,y;
        cin>>x>>y;
        t[x].push_back(y);
        t[y].push_back(x);
    }
    dep[s]=1;
    dfs(s,s);
    while(m--){
        int a,b;
        cin>>a>>b;
        cout<<lca(a,b)<<'\n';
    }
    return 0;
}
```
 - dsu on tree
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
 - Tarjan
```cpp
const int N=2e4+10;
vector<int> e[N];
int n,m,dfn[N],low[N],cnt,root,ans=0;
bool cut[N];
inline void tarjan(int u){
	dfn[u]=low[u]=++cnt;
	int rt_cnt=0;
	for(int v:e[u]){
		if(!dfn[v]){
			rt_cnt++;
			tarjan(v);
			low[u]=min(low[u],low[v]);
			if(low[v]>=dfn[u]&&u!=root) ans+=!cut[u],cut[u]=1;
		}
		else low[u]=min(low[u],dfn[v]);
	}
	if(u==root&&rt_cnt>1) ans+=!cut[u],cut[u]=1;
}
int main(){
	read(n),read(m);
	while(m--){
		int x,y;
		read(x),read(y);
		e[x].push_back(y);
		e[y].push_back(x);
	}
	for(int i=1;i<=n;i++) if(!dfn[i]) root=i,tarjan(i);
	write_enter(ans);
	for(int i=1;i<=n;i++) if(cut[i]) write_space(i);
	return 0;
}
```
 - 平板电视
```cpp
以下内容摘自lbw博客
Noip 开放了 pbds 的使用，这个东东可以支持很多奇妙的操作！
它可以支持哈希表，平衡二叉树，字典树（trie就行），堆（优先队列，功能非常之强大！
•	声明:
复制
// noi linux 中：
#include<bits/extc++.h>
using namespace __gnu_pbds;

// dev c++ 中：
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

________________________________________
它可以支持三种树，但似乎只有红黑树 rb_tree 不容易被卡，这里简单介绍一下它的一部分函数
•	insert(x) 在树中插入一个元素 x，返回 pair<iterator,bool>
•	erase(x) 在树种删除一个元素 x，返回 bool，表示删除操作是否成功
•	order_of_key(x) 返回 x 的排名
•	find_by_order(x) 返回排名为x的对应元素的迭代器（理解做指针
•	lower_bound(x)/upper_bound(x) 返回 定义序下的 第一个大于等于x/第一个大于'x' 的元素的迭代器
以下是其他的奇妙操作
•	join(x) ，将树x并入当前树，并把x删除
•	split(x,T2) 将树拆开，当前定义的序下小于等于x留在当前树，其余进入新树 T2
•	empty() 返回 bool 表示当前树是否为空
•	size() 返回整数表示树的大小
本题代码（左移是为了去重）：
复制
#include<bits/stdc++.h>
//#include<bits/extc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
#define ll long long
tree<ll, null_type, less<ll>, rb_tree_tag, tree_order_statistics_node_update> T;
int n;
ll k;
int main(){
	read(n);
	for(int i=1,opt;i<=n;i++){
		read(opt),read(k);
		if(opt==1) T.insert((k<<20)+i); //插入x
		if(opt==2) T.erase(T.lower_bound(k<<20)); //删除x
		if(opt==3) enter(T.order_of_key(k<<20)+1); //查询x的最小排名
		if(opt==4) enter(*T.find_by_order(k-1)>>20); //查询排名为x的数
		if(opt==5) enter(*--T.lower_bound(k<<20)>>20); //查询x的前驱
		if(opt==6) enter(*T.lower_bound((k+1)<<20)>>20); //查询x的后继
	}
	return 0;
}

本题另一种去重方法：
复制
#include<bits/stdc++.h>
//#include<bits/extc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
#define ll long long
tree<pair<ll,int>, null_type, less<pair<ll,int>>, rb_tree_tag, tree_order_statistics_node_update> T;
int n,cnt=0;
ll x;
int main(){
	read(n);
	for(int i=1,opt;i<=n;i++){
		read(opt),read(x);
		if(opt==1) T.insert({x,cnt++});
		if(opt==2) T.erase(T.lower_bound({x,0}));
		if(opt==3) enter(T.order_of_key({x,0})+1);
		if(opt==4) enter((*T.find_by_order(x-1)).first);
		if(opt==5) enter((*--T.lower_bound({x,0})).first);
		if(opt==6) enter((*T.lower_bound({x,cnt})).first);
	}
	return 0;
}
//idea offered by:cym
```
 - 整体二分
```cpp
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
```
 - 欧拉回路
所有点的度均为偶数的前提下，寻找欧拉回路，并把简单环拆分出来。

武汉邀请赛E
```cpp
void dfs(int u){
	for(int &i = cur[u]; i < t[u].size();){
		int v = t[u][i];
		if(e[HASH(u,v)] < E[HASH(u,v)]){
			e[HASH(u, v)]++, e[HASH(v, u)]++;
			i++; d[u]--, d[v]--;
			dfs(v);
		}
		else i++;
	}
	
	if(in_stk[u] == false){
		stk.push(u), in_stk[u] = true;
	}
	else{
		int pre_pos = u;
		col ++;
		while(!stk.empty()){
			int pos = stk.top();
			stk.pop();
			
			edge[HASH(pos, pre_pos)].push(col) , edge[HASH(pre_pos, pos)].push(col);
			pre_pos = pos;
			in_stk[pos] = false;
			if(pos == u) break;
		}
		
		stk.push(u), in_stk[u] = true;
	}
}
```
 - 莫比乌斯反演
 # 莫比乌斯反演

好久之前学习的内容了，但是似乎一直没有完全掌握，退役之前写一个博客。

 ## 零：需要知道的概念

  - 莫比乌斯函数

$$ \mu(n) =
\begin{cases}
1, \quad n=1 \\
0, \quad n含有平方因子 \\
(-1)^k, \quad n恰有k个不同的质因子
\end{cases}
$$

  - 莫比乌斯函数的性质 1
    
$$ \sum_{d|n} \mu(d) =
\begin{cases}
1, \quad n=1 \\
0, \quad n \neq 1
\end{cases}
$$ 

那么进一步地，一定有:
 $[gcd(i,j)==1] = \sum_{d|gcd(i,j)} \mu(d)$

  -  莫比乌斯函数的性质 2

$$ n = \sum_{d|n} \varphi(d) $$

 - 莫比乌斯函数的性质们

   莫比乌斯函数的性质 $1$ 可以写作
   
$$1 * \mu = \varepsilon$$

   莫比乌斯函数的性质 $2$ 可以写作
   
$$1 * \varphi = id$$

   类似的性质还有：

$$id^k * 1 = \sigma ^ k $$
$$ \mu * id = \varphi $$

  - 线性筛求莫比乌斯函数
```cpp
int mu[MAXN];
bool isnp[MAXN];
vector<int> primes;
void init(int n)
{
    mu[1] = 1;
    for (int i = 2; i <= n; i++)
    {
        if (!isnp[i])
            primes.push_back(i), mu[i] = -1; // 质数为-1
        for (int p : primes)
        {
            if (p * i > n)
                break;
            isnp[p * i] = 1;
            if (i % p == 0)
            {
                mu[p * i] = 0; // 有平方因数为0
                break;
            }
            else
                mu[p * i] = mu[p] * mu[i]; // 互质，利用积性函数性质
        }
    }
}
```

## 一：莫比乌斯变换

1. 形式一（常见形式）

$$已经知道 ： 1* \mu = \varepsilon$$

$$倘若有： f = 1 * g （那么总等价于 f * \varepsilon = 1 * g）$$

$$这时我们莫比乌斯反演得到： f * \mu = g （或者说： f * \mu = \varepsilon * g）$$

2. 形式二（神秘形式）

$$如果有 ： f(n) = \sum_{n|d} g(d)$$

$$那么一定有： g(n) = \sum_{n|d} f(d) \mu(\frac{d}{n})$$

3. 形式三（更加神秘的非卷积形式）
   
   对于数论函数 $f,g$ 和完全积性函数 $t$ 并且 $t(1)=1$ 而言:
   
$$
f(n) = \sum_{i=1}^{n} t(i) g\left(\left\lfloor \frac{n}{i} \right\rfloor\right)
$$

$$等价于$$

$$
\begin{align*}
g(n) = \sum_{i=1}^n \mu(i) t(i) f\left(\left\lfloor \frac{n}{i} \right\rfloor\right)
\end{align*}
$$

 - 自适应辛普森积分
```cpp
const db eps = 1e-10;
db a;
db F(db x){
	return pow(x, a/x-x);
};
db simpson(db a, db b){
	// calc [a,b]
	db c = (a + b)/2.0;
	return (b-a) * (F(a) + F(b) + F(c)*4.0 ) / 6.0;
}
db asr(db l, db r, db ans, int stp){
	if(stp < 0) return ans;
	db mid = (l+r)/2.0;
	db lft = simpson(l, mid), rgt = simpson(mid, r);
	if(fabs(lft+rgt - ans) < eps) return lft+rgt;
	else return asr(l, mid, lft, stp-1) + asr(mid, r, rgt, stp-1);
}

int main()
{
	cin>>a;
	if(a<0) cout<<"orz\n";
	else
		cout<<fixed<<setprecision(5)<<asr(eps, 20.0, simpson(eps, 20.0), 20)<<'\n';
	return 0;
}
```
 - Graham凸包
```cpp
template<class T>
class point {
	public:	
	T x, y, a;
	
	point() : x(T()), y(T()) {}
	point(T X, T Y) : x(X), y(Y) {}
	
	
	T operator ^ (const point &o) const{
		return x*o.y - y*o.x;
	}
	T operator * (const point &o) const{
		return x*o.x + y*o.y;
	}
	point operator * (const T &o) const{
		return point{x*o, y*o};
	}
	point operator + (const point &o) const{
		return (point<T>){x+o.x, y+o.y};
	}
	point operator - (const point &o) const{
		return (point<T>){x-o.x, y-o.y};
	}
	bool operator == (const point &o) const{
		return x == o.x && y == o.y;
	}
	
	T norm2() const { return x * x + y * y; }
	double norm() const { return hypot(x, y); }
	point rotate90() const { return {-y, x}; }
	point rotate(double ang) const {
		return {x * cos(ang) - y * sin(ang),
			x * sin(ang) + y * cos(ang)};
	}
	static bool collinear(const point &a, const point &b, const point &c) {
		return fabs((b - a) ^ (c - a)) < eps;
	}
};

template <class T>
class edge{
public:
	point<T> a, b, l;
	edge() : a(point<T>()), b(point<T>()), l(point<T>()) {}
	edge(point<T> a, point<T> b) : a(a), b(b){
		l = b - a; // a->b
	}
	bool operator == (const edge &o) const{
		return fabs(l ^ o.l) < eps;
	}
	T get_y (T x){
		if(fabs(a.x - b.x) < eps) return a.y;
		T res = a.y + (a.y-b.y)/(a.x-b.x) * (x - a.x);
		return res;
	}
	void print(string s = ""){
		cerr<<"edge["<<s<<"] : "<<a.x<<','<<a.y<<" -> "<<b.x<<','<<b.y<<'\n';
	}
};



void fst_IO()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);
}

int main()
{
	fst_IO();
	
	int n;
	cin>>n;
	vector< point<db> > p(n);
	for(int i=0;i<n;i++) cin>>p[i].x>>p[i].y;
	
	sort(p.begin(), p.end(), [&](point<db> a, point<db> b){
		if(a.y == b.y) return a.x < b.x;
		return a.y < b.y;
	});
	point<db> o = p[0];
//	cerr<<" point : "<<o.x<<' '<<o.y<<'\n';
	for(point<db> &cur : p){
		cur.a = atan2l(cur.y - o.y, cur.x - o.x);
		cur.y -= o.y, cur.x -= o.x;
	}
	sort(p.begin(), p.end(), [&](point<db> a, point<db> b){
		if(a.a == b.a) return a.norm2() < b.norm2();
		return a.a < b.a;
	});
	
	
	deque<int> stk;
	stk.push_back(0);
	for(int i=1;i<n;i++){
		while( stk.size()>=2 ){
			point<db> p1 = p[stk[stk.size()-1]], p2 = p[stk[stk.size()-2]];
			point<db> vec1 = p[i] - p1, vec2 = p1 - p2;
			if (p1 == p[i] || (vec1 ^ vec2) >= 0) stk.pop_back();
			else break;
		}
		stk.push_back(i);
	}
	
	db ans = 0;
	for(int i=0;i<stk.size();i++){
		point<db> vec = p[stk[(i+1)%stk.size()]] - p[stk[i]];
		ans += vec.norm();
	}
	
	cout<<fixed<<setprecision(2)<<ans<<'\n';
	return 0;
}
```
 - 半平面交
```cpp
const db eps = 1e-7;

template<class T>
class point {
	public:	
	T x, y;
	
	point() : x(T()), y(T()) {}
	point(T X, T Y) : x(X), y(Y) {}
	
	
	T operator ^ (const point &o) const{
		return x*o.y - y*o.x;
	}
	T operator * (const point &o) const{
		return x*o.x + y*o.y;
	}
	point operator * (const T &o) const{
		return point{x*o, y*o};
	}
	point operator + (const point &o) const{
		return (point<T>){x+o.x, y+o.y};
	}
	point operator - (const point &o) const{
		return (point<T>){x-o.x, y-o.y};
	}
	bool operator == (const point &o) const{
		return x == o.x && y == o.y;
	}
	
	// 以下：长度、长度平方、旋转、判共线
	T norm2() const { return x * x + y * y; }          // 长度平方
	double norm() const { return hypot(x, y); }        // 长度
	point rotate90() const { return {-y, x}; }         // 逆时针 90°
	point rotate(double ang) const {                   // 任意角（弧度）
		return {x * cos(ang) - y * sin(ang),
			x * sin(ang) + y * cos(ang)};
	}
	static bool collinear(const point &a, const point &b, const point &c) {
		return fabs((b - a) ^ (c - a)) < eps;
	}
};

template<class T>
point<T> intersection_point(const point<T> &s1,const point<T> &t1,const point<T> &s2,const point<T> &t2){
	point<T> l1=s2-s1, l2=t2-s1, l3=s2-t1, l4=t2-t1;
	db tmp1=l1^l2, tmp2=l3^l4;
	db k = (tmp1)/(tmp1-tmp2);
	point<T> l0 = t1-s1;
	l0 = l0*k;
	point<T> res=s1;
	res.x+=l0.x, res.y+=l0.y;
	return res;
}

template <class T>
class edge{
public:
	point<T> a, b, l;
	edge() : a(point<T>()), b(point<T>()), l(point<T>()) {}
	edge(point<T> a, point<T> b) : a(a), b(b){
		l = b - a; // a->b
	}
	bool operator == (const edge &o) const{
		return fabs(l ^ o.l) < eps;
	}
	T get_y (T x){
		if(fabs(a.x - b.x) < eps) return a.y;
		T res = a.y + (a.y-b.y)/(a.x-b.x) * (x - a.x);
		return res;
	}
	void print(string s = ""){
		cerr<<"edge["<<s<<"] : "<<a.x<<','<<a.y<<" -> "<<b.x<<','<<b.y<<'\n';
	}
};

template <class T>
class convex{
public:
	vector<edge<T>> e;
	convex() {}
	convex(vector<edge<T>> e) : e(e) {}
	
};

template<class T>
vector<edge<T>> HalfPlane_intersection(vector<edge<T>> e, db minx, db maxx){
	sort(e.begin(), e.end(), [&](edge<T> e1, edge<T> e2){
		db t1 = atan2(e1.l.y, e1.l.x);
		db t2 = atan2(e2.l.y, e2.l.x);
		
		if(fabs(t1 - t2) > eps) return t1 < t2;
		else{	// 应该不能简单判断横坐标大小来判断
			point<T> e3 = e1.a-e1.a, e4 = e2.b-e1.a;
			return (e3 ^ e4) > 0;  // > eps ? 
		}
	});
	
	
	auto check = [](const edge<T> &a, const edge<T> &b, const edge<T> &c) -> bool {
		// 判断b,c交点是否在直线a右侧
		point<T> p = intersection_point(b.a, b.b, c.a, c.b);
		point<T> tmp = p - a.a;
		return ( tmp ^ a.l ) >= 0;
	};
	
	int n = e.size();
	vector<edge<T>> q(n);
	q[0] = e[0];
	int head = 0, tail = 0;
	for(int i=1;i<n;i++){
		if(e[i] == e[i-1]) continue;
		if(head < tail){
			point<T> p = intersection_point(q[tail].a, q[tail].b, q[tail-1].a, q[tail-1].b);
		}
		while(head < tail && check(e[i], q[tail], q[tail-1])){
			tail--;
		}
		while(head < tail && check(e[i], q[head], q[head+1])){
			head++;
		}
		q[++tail] = e[i];
	}
	
	while(head < tail-1 && check(q[head],q[tail],q[tail-1])) tail--;
	vector<edge<T>> result;
	for(int i=head;i<=tail;i++){
		point<T> p1(minx, q[i].get_y(minx)), p2(maxx, q[i].get_y(maxx));
		q[i] = edge<T>(p1, p2);
	}
	for(int i=head;i<=tail;i++){
		if(i+1<=tail){
			point<T> p = intersection_point(q[i].a, q[i].b, q[i+1].a, q[i+1].b);
			q[i].b = q[i+1].a = p;
		}
		if(fabs(q[i].a.x - q[i].b.x) < eps) continue;
		result.push_back(q[i]);
	}
	
	return result;
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
	int n;
	cin>>n;
	vector<point<db>> p(n);
	
	for(int i=0;i<n;i++) cin>>p[i].x;
	for(int i=0;i<n;i++) cin>>p[i].y;
	
	db ans = 1e21;
	for(int splt = -1; splt <= n; splt ++){
		vector<edge<db>> e;
		for(int i=0;i+1<splt;i++){
			e.push_back(edge<db>(p[i] , p[i+1]));
		}
		for(int i=max(0, splt);i+1<n;i++){
			e.push_back(edge<db>(p[i], p[i+1] ));
		}
		
		e = HalfPlane_intersection(e, p[0].x, p[n-1].x);
		vector<db> X;
		auto check_leagal = [&](db x){
			if(splt == -1) return x <= p[0].x;
			if(splt == n)  return x >= p[n-1].x;
			return x <= p[splt].x && x >= p[splt-1].x;
		};
		for(int i=0;i<n;i++) if(check_leagal(p[i].x)) X.push_back(p[i].x);
		for(edge l : e){
			if(check_leagal(l.a.x)) X.push_back(l.a.x);
			if(check_leagal(l.b.x)) X.push_back(l.b.x);
		}
		
		for(db x: X){
			db y = -1;
			for(int i=0;i+1<n;i++) if(p[i].x <= x && x<= p[i+1].x){
				edge<db> E(p[i], p[i+1]);
				y = E.get_y(x);
				break;
			}
			if(y < 0) continue;
			
			for(edge l:e) if(l.a.x <= x && x <= l.b.x){
				cmin(ans, l.get_y(x) - y);
			}
		}
	}
	cout<<fixed<<setprecision(3)<<ans<<'\n';
	return 0;
}
```
 - 计算几何板子合集
```cpp
#include<bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db double
#define ldb long double
#define max(A,B) (A>B? A:B)
#define min(A,B) (A<B? A:B)
#define yes cout<<"Yes\n"
#define no cout<<"No\n"
#define pir(X) pair<X,X>
#define mpr(A,B) make_pair(A,B)
#define fr first
#define sc second
#define sq(x) ((x)*(x))
#define sgn(x) (((x)>eps) - ((x)<-eps))
using namespace std;

const db eps = 1e-7;
const db pii=3.1415926;

class Point{
public:
    db x,y,arctan;
    Point(db _x=0,db _y=0) : x(_x),y(_y) {}
    void read() { cin>>x>>y; }
    void print(string s=""){
        cerr<<" Point ["<<s<<"] : x="<<x<<" y="<<y<<'\n';
    }

    Point turn(const db k) const {
        Point p(x*cos(k)-y*sin(k),x*sin(k)+y*cos(k));
        return p;
    }
    Point operator - (const Point &o) const{
        Point ret(x-o.x,y-o.y);
        return ret;
    }
    Point operator + (const Point &o) const{
        Point ret(x+o.x,y+o.y);
        return ret;
    }
    bool operator == (const Point &o) const{
        return x==o.x&&y==o.y;
    }
};
db dis(Point a,Point b){
    return sqrt(sq(a.x-b.x)+sq(a.y-b.y));
}

class Vector2{
public:
    db x,y;
    Vector2(db _x=0,db _y=0): x(_x),y(_y) {}
    void read() { cin>>x>>y; }
    void print(string s=""){
        cerr<<" Vector2 ["<<s<<"] : x="<<x<<" y="<<y<<'\n';
    }

    void assign(db _x,db _y) { x=_x, y=_y; }
    void assign(const Point &a,const Point &b){ // 向量赋值，注意前后形参的顺序
        x=a.x-b.x;
        y=a.y-b.y;
    }
    Vector2 operator + (const Vector2 &o) const {
        return (Vector2){x+o.x,y+o.y};
    }
    Vector2 operator - (const Vector2 &o) const {
        return (Vector2){x-o.x,y-o.y};
    }
    Vector2 operator * (const db a) const {
        return (Vector2){a*x,a*y};
    }
    db operator * (const Vector2 &o) const {
        return x*o.x+y*o.y;
    }
    db operator ^ (const Vector2 &o) const {
        return x*o.y-y*o.x;
    }
    db len() const{
        return sqrt(x*x+y*y);
    }
    db len2() const{
        return x*x+y*y;
    }
    Vector2 unit() const {
        db len=sqrt(x*x+y*y);
        Vector2 l(x/len,y/len);
        return l;
    }
};

Point Projection(const Point& p,const Point& p1,const Point& p2){ // 计算 p 点到 p1-p2 上的投影点坐标
    Vector2 l0(p.x-p1.x,p.y-p1.y);
    Vector2 l1(p2.x-p1.x,p2.y-p1.y);
    db tmp=l0*l1/l1.len();
    Vector2 l2=l1.unit()*tmp;
    Point ret(p1.x+l2.x,p1.y+l2.y);
    return ret;
}
Point Reflection(const Point& p,const Point& p1,const Point& p2){ // p 关于 p1-p2 的对称点
    Point p0=Projection(p,p1,p2);
    Point ret(p0.x*2-p.x,p0.y*2-p.y);
    return ret;
}
bool Parallel(const Vector2 &l1,const Vector2 &l2){ // 判断平行
    return fabs(l1^l2)<=eps;
}
bool Orthogonal(const Vector2 &l1,Vector2 &l2){ // 判断垂直/正交
    return fabs(l1*l2)<=eps;
}
bool In_line(const Point &a,const Point &b,const Point &c){  // 判断 b 是否在线段 a-c 上
    if(!(min(a.x,c.x)<=b.x+eps&&b.x<=max(a.x,c.x)+eps)) return false;
    if(!(min(a.y,c.y)<=b.y+eps&&b.y<=max(a.y,c.y)+eps)) return false;
    Vector2 l1,l2;
    l1.assign(a,b),l2.assign(c,b);
    return fabs(l1^l2)<=eps;
};
namespace Intersection{

    /*
        判断两条线段是否有交点
        求两条直线交点
    */

    bool interval_intersect_check(db l1,db r1,db l2,db r2){  // 判断 [l1,r1] 和 [l2,r2] 两个区间是否有交
        if(l1>r1) swap(l1,r1);
        if(l2>r2) swap(l2,r2);
        if(l1<=l2+eps&&l2<=r1+eps) return true;
        if(l2<=l1+eps&&l1<=r2+eps) return true;
        return false;
    };
    bool divided_by_line(const Point &a,const Point &b,const Point &p1,const Point &p2){ // 判断p1-p2是否被a-b分割(对于p2而言)
        Vector2 l0,l1,l2;
        l0.assign(p1,p2),l1.assign(a,p2),l2.assign(b,p2);
        db tmp1=l0^l1,tmp2=l0^l2;
        return (tmp1<0&&tmp2>0)||(tmp1>0&&tmp2<0);
    };

    bool Intersection_judge(const Point &s1,const Point &t1,const Point &s2,const Point &t2){ // 判断s1-t1 , s2-t2两条线段是否有交
        if(In_line(s1,s2,t1)||In_line(s1,t2,t1)||In_line(s2,s1,t2)||In_line(s2,t1,t2)) return 1;
        Vector2 a,b;
        a.assign(t1,s1),b.assign(t2,s2);
        if(Parallel(a,b)) return 0;
        if(interval_intersect_check(s1.x,t1.x,s2.x,t2.x)&&interval_intersect_check(s1.y,t1.y,s2.y,t2.y)){
            if(!divided_by_line(s1,t1,s2,t2)) return 0;
            // if(!divided_by_line(s1,t1,t2,s2)) return 0;
            if(!divided_by_line(s2,t2,s1,t1)) return 0;
            // if(!divided_by_line(s2,t2,t1,s1)) return 0;

            return 1;
        }
        else return 0;
    }
    Point Intersection_point(const Point &s1,const Point &t1,const Point &s2,const Point &t2){
        Vector2 l1,l2,l3,l4;
        l1.assign(s2,s1),l2.assign(t2,s1);
        l3.assign(s2,t1),l4.assign(t2,t1);
        db tmp1=l1^l2,tmp2=l3^l4;
        db k=(tmp1)/(tmp1-tmp2);
        Vector2 l0;l0.assign(t1,s1);
        l0=l0*k;
        Point res=s1;
        res.x+=l0.x,res.y+=l0.y;
        return res;
    }
}

db dot_to_segment(const Point &p,const Point &a,const Point &b){    // p 到线段 a-b 的距离
    Point pp=Projection(p,a,b);
    Vector2 l1,l2;
    l1.assign(a,p),l2.assign(b,p);
    if(In_line(a,pp,b)){
        Vector2 l;
        l.assign(b,a);
        return fabs(l1^l2)/(l.len());
    }
    else return min(l1.len(),l2.len());
}
db dot_to_line(const Point &p,const Point &a,const Point &b){    // p 到直线 a-b 的距离
    Vector2 l1,l2;
    l1.assign(a,p),l2.assign(b,p);
    Vector2 l;
    l.assign(b,a);
    return fabs(l1^l2)/(l.len());
}


db Distance(const Point &a1,const Point &b1,const Point &a2,const Point &b2){   // 线段a1-b1到线段a2-b2的距离
    if(Intersection::Intersection_judge(a1,b1,a2,b2)) return 0;
    return min( min(dot_to_segment(a1,a2,b2),dot_to_segment(b1,a2,b2)) , min(dot_to_segment(a2,a1,b1),dot_to_segment(b2,a1,b1)) );
}

bool comp1(Point a,Point b){    // Andrew 求凸包用到的排序
    if(a.x==b.x) return a.y>b.y;
    else return a.x>b.x;
}

class Polygon{
public:
    int n;
    vector<Point> pt;
    Polygon(){ n=0; pt.clear(); }
    Polygon(int siz) :n(siz) { pt.resize(siz); }
    Polygon(const vector<Point> &p){
        n=p.size();
        pt=p;   // 大概是这样？
        // pt.resize(n);
        // for(int i=0;i<n;i++) pt[i]=p[i];
    }
    void read(){
        cin>>n;
        pt.resize(n);
        for(int i=0;i<n;i++) pt[i].read();
    }
    void print(string s=""){
        cerr<<" Polygon ["<<s<<"] : n="<<n<<'\n';
        for(Point p:pt) cerr<<"     x="<<p.x<<" y="<<p.y<<'\n';
    }

    db Square(){    // 求面积(不需要形状是凸的) 前提是 pt 中的点是逆时针顺序
        Point o(0,0);
        db S=0;
        for(int i=0;i<n;i++){
            Vector2 l1,l2;
            l1.assign(pt[i],o),l2.assign(pt[(i+1)%n],o);
            S+=l1^l2;
        }
        return S/2;
    }
    db Perimeter(){     // 求凸包周长(前提是pt中的点是有序给出的)
        db C=0;
        for(int i=0;i<n;i++){
            Vector2 l1;
            l1.assign(pt[i],pt[(i+1)%n]);;
            C+=l1.len();
        }
        return C;
    }
    bool Convex_judge(){  // 判断是否是凸多边形 前提是 pt 中的点是逆时针顺序
        for(int i=0;i<n;i++){
            Vector2 l1,l2;
            l1.assign(pt[(i+1)%n],pt[i]),l2.assign(pt[(i+2)%n],pt[(i+1)%n]);
            if((l1^l2)<0) return 0;
        }
        return 1;
    }
    bool intersect_judge2(const Point &a,Point b,Point c){  // 判断 y=a.y 是否与线段 b-c(左闭右开) 有交，这里不写 '&'
        if(b.y<c.y) swap(b,c);
        if(!(c.y<=a.y+eps&&a.y<b.y)) return false;
        db x=b.x+(b.y-a.y)*(c.x-b.x)/(b.y-c.y);
        return x+eps>=a.x;
    }
    int In_Convex(Point p){ // O(n)判断点是否在凸包(多边形？)内/  0:在外部 ;  1:在边上 ; 2:在内部
        for(int i=0;i<n;i++) if(In_line(pt[i],p,pt[(i+1)%n])) return 1;
        int cnt=0;
        for(int i=0;i<n;i++) cnt+=(int)intersect_judge2(p,pt[i],pt[(i+1)%n]);
        if(cnt&1) return 2;
        else return 0;
    }
    
    Polygon Convex_Hull(){ // Andrew 法求(最大)凸包，这里会改变 pt 内顺序
        Point oo=pt[0];
        for(int i=0;i<n;i++) if(pt[i].y<oo.y||(pt[i].y==oo.y&&pt[i].x<oo.x)) oo=pt[i];
        sort(pt.begin(),pt.end(),comp1);
        int tot=0;
        vector<Point> stk;
        stk.resize(n+2);
        for(int i=0;i<n;i++){
            while(tot>1){
                Vector2 l1,l2;
                l1.assign(pt[i],stk[tot-2]),l2.assign(stk[tot-1],stk[tot-2]);
                if(sgn(l1^l2)>=0) tot--;    // 这里可以加一个等号，避免重复点，在卡壳中有用(特判掉只有一个点等情况)
                else break;
            }
            stk[tot++]=pt[i];   // 注意这里是 tot++
        }
        int lim=tot;
        for(int i=n-2;i>=0;i--){
            while(tot>lim){
                Vector2 l1,l2;
                l1.assign(pt[i],stk[tot-2]),l2.assign(stk[tot-1],stk[tot-2]);
                if(sgn(l1^l2)>=0) tot--;
                else break;
            }
            stk[tot++]=pt[i];
        }
        if(n>1) tot--;
        while(stk.size()>tot) stk.pop_back();
        Polygon res;
        res.n=tot,res.pt.resize(tot);
        
        /*
            这里的凸包会满足一些性质：
                边数最多，pt[0]一定是y最小且y相等时x最小的点，pt内的点按照逆时针顺序给出
        */
        int id;
        for(int i=0;i<tot;i++) if(stk[i].x==oo.x&&stk[i].y==oo.y) id=i;
        for(int i=0;i<tot;i++) res.pt[i]=stk[(i+id)%tot];
        return res;
    }
    
    db Roating_Calipers(){  // 求凸包直径 (前提是已经是凸包)
        db ans=0;
        if(n==1) return 0;
        if(n==2) return dis(pt[0],pt[1]);
        vector<int> stk(n+2);
        for(int i=0;i<n;i++) stk[i]=i;
        int top=n;
        stk[++top]=0;
        int cur=0;
        for(int i=0;i+1<top;i++){
            while(1){
                int nxt=cur+1;
                if(nxt>=top) nxt-=top;
                if( dot_to_line(pt[stk[cur]],pt[stk[i]],pt[stk[i+1]]) <= dot_to_line(pt[stk[nxt]],pt[stk[i]],pt[stk[i+1]]) + eps ) cur=nxt;
                else break;
            }
            ans=max(ans, max(dis(pt[stk[i]],pt[stk[cur]]), dis(pt[stk[i+1]],pt[stk[cur]]) ) );
        }
        
        return ans;
    }

    int check_k(const Point &p,const Point &p1,const Point &p2){  // 判断： p-p2的斜率大于p-p1的斜率
        Vector2 tmp1,tmp2;
        tmp1.assign(p1,p);
        tmp2.assign(p2,p);
        return sgn(tmp1^tmp2);
    }
    bool In_Convex2(const Point &p){    // O(logn) 时间判断点是否在凸包内,要求凸包上的点逆时针给出,pt[0]是y最小的点
        if(check_k(pt[0],p,pt[1])==1) return 0;
        if(check_k(pt[0],p,pt[1])==0) return In_line(pt[0],p,pt[1]);
        if(check_k(pt[0],pt[n-1],p)==1) return 0;
        if(check_k(pt[0],pt[n-1],p)==0) return In_line(pt[0],p,pt[n-1]);
        int l=1,r=n-2,rec,mid;
        while(l<=r){
            mid=(l+r)>>1;
            if(check_k(pt[0],pt[mid],p)>=0) rec=mid,l=mid+1;
            else r=mid-1;
        }
        int res=check_k(p,pt[rec],pt[rec+1]);
        return res>=0;
    }
    vector<Vector2> convex_vector() const {    // 生成一组凸包上的边，用于闵可夫斯基求和
        vector<Vector2> ret(n);
        for(int i=0;i<n;i++) ret[i].assign(pt[(i+1==n? 0:i+1)],pt[i]);
        return ret;
    }
    
};

Polygon Minkowski (const Polygon &o1,const Polygon &o2){   // 闵可夫斯基和
    vector<Vector2> s1=o1.convex_vector(),s2=o2.convex_vector();
    vector<Point> p(o1.n+o2.n+1);
    int tot=0,l1=0,l2=0;
    p[0]=o1.pt[0]+o2.pt[0];
    while(l1<o1.n&&l2<o2.n){
        Point tmp=p[tot++];
        if((s1[l1]^s2[l2])>0) tmp.x+=s1[l1].x,tmp.y+=s1[l1].y,l1++;
        else tmp.x+=s2[l2].x,tmp.y+=s2[l2].y,l2++;
        p[tot]=tmp;
    }
    while(l1<o1.n){
        Point tmp=p[tot++];
        tmp.x+=s1[l1].x,tmp.y+=s1[l1].y,l1++;
        p[tot]=tmp;
    }
    while(l2<o2.n){
        Point tmp=p[tot++];
        tmp.x+=s2[l2].x,tmp.y+=s2[l2].y,l2++;
        p[tot]=tmp;
    }
    Polygon ret(p);
    ret=ret.Convex_Hull();  // 重做一遍凸包，防止奇怪的问题
    return ret;
}

void fst_IO(){
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);
}

int main(){
    fst_IO();
    int n,m,q;
    cin>>n>>m>>q;
    Polygon o1(n),o2(m);
    for(int i=0;i<n;i++) o1.pt[i].read();
    Point oo(0,0);
    for(int i=0;i<m;i++) o2.pt[i].read(),o2.pt[i]=oo-o2.pt[i];
    o1=o1.Convex_Hull();
    o2=o2.Convex_Hull();
    Polygon poly=Minkowski(o1,o2);
    while(q--){
        Point p;
        p.read();
        cout<<poly.In_Convex2(p)<<'\n';
    }
    return 0;
}
```
