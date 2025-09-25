[P2600 [ZJOI2008] 瞭望塔](https://www.luogu.com.cn/problem/P2600)

求开放半平面的半平面交

```
#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define db long double
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
