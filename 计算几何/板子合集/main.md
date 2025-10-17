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
