```cpp
#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db double
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

luogu : [https://www.luogu.com.cn/problem/P4526](https://www.luogu.com.cn/problem/P4526)
# P4526 【模板】自适应辛普森法 2

试计算积分

$\displaystyle{\int_0^\infty x^{\frac{a}{x}-x}\mathrm{d}x}$

保留至小数点后 $5$ 位。若积分发散，请输出 $\text{orz}$ 。
