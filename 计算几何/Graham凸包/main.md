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
