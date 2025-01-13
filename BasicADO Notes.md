## BasicADO Notes

图 $G=(V,E)$，$E=\{(u,v,w)\}$，$n=|V|,m=|E|$.

#### 采样

取 $k=3$，初始化 $A_{0}=V$.

对于 $A_i$，每个点 $x\in A_{i-1}$ 有 $n^{-\frac 1 k}$ 采样成功加入 $A_i$，最终 $A_i$ 的期望大小为 $n^{1-\frac i k}$，直接钦定 $A_k=\empty$，其他的 $A$ 采样出空重新采样即可。

#### 求 C/B

约定数组：$B(x)=\{u,d\}$，表示最终需要求得的 DO，需要 $(x,u)\rightarrow d$ 的随机索引，或者返回不存在。

$k\times n$ 的数组  $P(i,x)=(rt,dis)$ 表示 $A_i$ 为起点集合得到的最短路树，表示 $x$ 到 $A_i$ 的最近点 $rt$，距离为 $dis$，约定在距离相等时首先选择在 $j$ 更大的 $A_{j}$ 中的点，其次约定最小编号。

从 $k-1$ 倒序到 $0$ 逐步进行如下过程计算：

1. 求 $A_i$ 点为源点集合的最短路树，得到 $P(i,\dots)$ 数组。 使用 pq 实现，本步复杂度为 $O(m\log n)$。

2. 以 $A_{i+1}-A_i$ 中的点 $x$ 为源点分别跑一次局部的 dijk，求数组 $C[x]$ 表示 $x$​ 点为根的最短路树信息。

   - 关于局部到何时进行剪枝？

     如果 $dis(x\rightarrow y) \ge P(i+1,y)$ 就直接减掉子树。

3. 当前层 $i$ 求得的 $C[rt]=\{x,d\}$ 是 $B(x)=\{rt,d\}$ 的一部分与其他地方不重复的反数组，将其加入 $B$。

论文证明了所有的 $C$ 经过剪枝后没有重复对 $B$ 产生贡献，且 $B$ 的期望大小是 $O(kn^{1+\frac 1 k})$.

因此 2 3 部分的实现作为时间复杂度瓶颈，应为 $O(kn^{1+\frac 1 k}\log n)$，但是这里用 pq 可能达不到。

#### 查询

要求连通。

```cpp
DistanceType Query(NodeType u, NodeType v)
{
    int i = 0;
    NodeType w = u;
    while(B[v].count(w) == 0)
    {
        i++;
        std::swap(u, v);
        w = P[i][u].first;
    }
    return B[u][w] + B[v][w];
}
```





