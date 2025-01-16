#### 总述

本文是关于 [2005]Approximate Distance Oracles 的一篇实践性中文介绍。

解决的问题是经过对原图的预处理，存储一定空间大小的结构（称 Distance Oracles、DO），然后快速回答原图任意两点（近似）最短路的结构。特别的，要求图为无向带非负权连通图。

如其名，文中给出了一种空间复杂度 $O(kn^{1+\frac 1 k})$ 的 DO ，构建复杂度大致为结构复杂度乘 $\log$，回答的两点距离具有 $2k-1$ 近似比（即实际查询结果在真实结果 $2k-1$ 倍内），可以改进为路径输出的模式而不是仅回答距离。这里的 $k$ 为大于 $1$​ 的整数。 

考虑实际算法卡常时，主要关注结构的空间复杂度和实际近似比。

#### 约定

输入 $G=(V,E)$，$E=\{(u,v,w)\}$，$n=|V|,m=|E|$​​.

输出 Label $B(v)=\{u,d\}$，即对于图中每一个点维护 $v$ 维护一个集合，表示 $v$ 到 $u$ 的最短路径为 $d$，结构需要支持回答 $u$ 是否在 $B(v)$ 中；辅助结构 $P[i][x]$ 为一个 $k\times n$ 的数组，具体内容稍后介绍。

#### 构建过程

在决定 $k$ 后，首先采样出集合 $A_0,A_1,A_2,\dots,A_k$

令 $A_0=V$，然后对于 $A_i$，每个点 $x\in A_{i-1}$ 有 $n^{-\frac 1 k}$ 采样成功加入 $A_i$，最终 $A_i$ 的期望大小为 $n^{1-\frac i k}$，直接钦定 $A_k=\empty$。其他的 $A$ 若采样出空重新采样直到非空即可。

随机采样的过程在后续证明提供了两点作用，其一是期望大小，其二是 $\mathtt {Pr}[w\in A_{i+1}|w\in A_i]=n^{-\frac 1 k}$.

接下来从 $k-1$ 到 $0$ 倒序求解结构：

1. 求 $A_i$ 中所有点为源点集合的最短路树，得到 $P[i][x]$ 数组，表示 $x$ 距离 $A_i$ 中哪个节点最近以及具体距离。

2. 以 $A_i-A_{i+1}$ 中的每个点 $u$ 为源点分别跑一次**局部的** dijk，求临时集合 $C(u)=\{v,d\}$ 表示 $u$​ 点为根的最短路树信息。

   局部：在 $u$ 为根的 dijk 中，如果当前距离 $dis(u\rightarrow y)\ge P[i+1][y]$，那么整颗子树都剪掉。

3. 枚举所有 $\{v,d\} \in C(u)$ ，将 $\{u,d\}$ 加入 $B(v)$。

#### 查询

假设需要查询 $u,v$。初始 $u_0=u,v_0=v$。

从 $i=0$ 开始，令 $w_i=p[i][u_i]$。

若 $w_i\in B(v_i)$ ，则返回 $v_i,w_i$ 在 $B(v_i)$ 的距离 $d_{v_i,w_i}$ 和类似的 $d_{u_i,w_i}$；

否则令 $i\leftarrow i+1,u_{i+1}\leftarrow v_i,v_{i+1}\leftarrow u_i$。

#### 证明内容

1. $\sum B(v)\sim O(kn^{1+\frac 1 k})$ 

   对于每个 $B(v)$，证明对每个 $0\le i \le  k-1$ 有 $\mathbb E[B(v)\cap A_i]\le n^{\frac 1 k}$，因为 $A_{k-1}\sim n^{\frac 1 k}$，因此对于 $i=k-1$ 显然成立。

   对于 $i<k-1$ 的情况，将 $w_1,w_2,\dots,w_l\in A_{i}$ 按距离 $v$ 的顺序排序，若某个 $w_j\in B(v)$，说明 $\delta(w_j,v)<\delta(A_{i+1},v)$ 那么显然 $w_1,w_2,\dots,w_{j-1}\notin A_{i+1}$。由于 $\mathtt {Pr}[w\in A_{i+1}|w\in A_i]=n^{-\frac 1 k}$，有 $\mathtt{Pr}[w_j\in B(v)]\le (1-n^{-\frac 1 k})^{j-1}\le n^{\frac 1 k}$​.

   这里说明 $B$ 期望大小与条件概率和 $A$ 的期望大小有关。

2. 查询会停止且近似比为 $2k-1$​.

   因为 $A_{k-1}\subseteq B(v)$，因此如果前面没有停止，$k-1$ 一定存在中间节点。

   注意 $w_i=p[i][u_i]$，因此可以正确查询距离 $\delta(w_i,u_i)$，记正确距离 $\Delta=\delta(u_0,v_0)$，接下来证明 $\delta(w_i,u_i)\le\delta(w_{i-1},u_{i-1})+\Delta$

   因为 $w_{i-1}\notin B(v_{i-1})$ 所以 $\delta(w_{i-1},v_{i-1})\ge \delta(A_i,v_{i-1})=\delta(p[i][u_i],v_{i-1})=\delta(w_i,u_i)$

   由三角形不等式 $\delta(w_i,u_i)\le \delta(w_{i-1},v_{i-1})\le \delta(w_{i-1},u_{i-1})+\Delta$，得证。

   因此 $\delta(w,u)\le (k-1)\Delta$，同样三角形不等式 $\delta(w,v)\le (k-1)\Delta+\Delta$，故近似比为 $2k-1$.

   这里说明近似比之和迭代层数有关，与选点是否随机无关。

   

   

   

   

   

   

