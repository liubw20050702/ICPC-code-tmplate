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

## 二：莫比乌斯反演的本质
