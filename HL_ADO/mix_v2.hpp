#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <limits>
#include <cmath>
#include <random>
#include "../AdjacencyList.hpp"

using std::vector;
using std::pair;
using std::sort;
using std::min;
using std::swap;

template<typename NodeType, typename DistanceType>
class Mix_v2
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    private:
        AdjList &graph;
        vector<vector<pair<NodeType, DistanceType>>> L;
        vector<DistanceType> qtable;
        DistanceType inf;
        int seed, K;

    public:
        long long LabelSize;
        void setseed(int Seed){seed=Seed;}
        Mix_v2(AdjList &lis, int K): graph(lis), K(K)
        {
            LabelSize = 0;
            std::cerr << "Start Mix_v1!" << std::endl;
            NodeType n=graph.GetSize();
            inf=std::numeric_limits<DistanceType>::max();

            L.assign(n, vector<pair <NodeType, DistanceType> >());
            qtable.assign(n, inf);

            vector <int> deg(n);
            for(NodeType u=0;u<n;u++)
                for(auto v: graph.GetAllEdges(u))
                    deg[v.first]++;

            vector <NodeType> vid(n);
            for(NodeType i=0;i<n;i++) vid[i]=i;
            sort(vid.begin(),vid.end(),[&](auto x,auto y){return deg[x]>deg[y];});

            vector <int> Siz(K + 1);
            for(int i = 0; i < K; i++) Siz[i] = ceil(pow(graph.GetSize(), 1.0 - 1.0*i/K));

            int sK = Siz[K-1];
            vector <int> A0;
            for(int i = sK; i < n; i++) A0.push_back(vid[i]);

            std::mt19937 gen(seed);
            shuffle(A0.begin(), A0.end(), gen);
            for(int i = sK; i < n; i++) vid[i] = A0[i - sK];

            vector <DistanceType> d(n,inf);

            //for(int i = 0; i < K; i++) std::cerr << Siz[i] <<" ";std::cerr<<"\n";
            
            int curk = K-1;
            int lasSize = 0;
            vector <DistanceType> pi(n,inf);
            vector <NodeType> isr(n);
            for(int i=0;i<=n;i++)
            {
                if(i == Siz[curk] || i == n)
                {
                    int psize = 0;
                    for(int j = 0; j < n; j++) psize += L[j].size();
                    std::cerr << "Curk = " << curk << ", vextices Label Size: " << psize - lasSize << std::endl;
                    lasSize = psize;

                    if(i==n) break;
                    
                    vector <NodeType> Ai;
                    for(int j = 0; j < i; j++) Ai.push_back(vid[j]);
                    pi = SSSPTree(Ai);
                    --curk;
                }
                vector <NodeType> inq;
                NodeType rt=vid[i];
                d[rt]=0;

                for(auto [v, w]: L[rt])
                    qtable[v] = w;
                int flg = 1;

                std::priority_queue <pair <DistanceType, NodeType>, vector <pair <DistanceType, NodeType> >, std::greater <pair <DistanceType, NodeType>> > q;
                q.push({d[rt],rt});
                while(!q.empty())
                {
                    auto [dis, u]=q.top();
                    q.pop();
                    if(d[u]<dis) continue;
                    inq.push_back(u);

                    if(isr[u]) 
                    {
                        flg = 0;
                        continue;
                    }

                    if(flg || Query_Greater(u, d[u]))
                    {
                        L[u].push_back({rt,d[u]});
                        for(auto [v,w]: graph.GetAllEdges(u))
                            if(d[u]+w<d[v] && d[u]+w<pi[v])
                            {
                                d[v]=d[u]+w;
                                q.push({d[v],v});
                            }
                    }
                }
                for(auto u: inq) d[u]=inf;
                for(auto [v, _]: L[rt])
                    qtable[v] = inf;
                isr[rt] = 1;
            }

            for(int i = 0; i < n; i++) LabelSize += L[i].size();
            std::cerr << "HL Contruct Done!" << std::endl;
        }
        
        bool Query_Greater(NodeType u, DistanceType d)
        {
            for(auto [v, w]: L[u]) if(qtable[v] <= d - w) return false;
            return true;
        }

        std::vector<DistanceType> SSSPTree(std::vector <int> &Ai)
        {
            std::vector<DistanceType> dist(graph.GetSize(), inf);

            std::priority_queue<std::pair<DistanceType, NodeType>, std::vector<std::pair<DistanceType, NodeType>>, std::greater<std::pair<DistanceType, NodeType>>> q;

            std::vector<std::pair<NodeType, DistanceType>> res(graph.GetSize(), {-1, inf});
            std::vector<DistanceType> ret(graph.GetSize(), inf);

            for(NodeType x: Ai)
            {
                dist[x] = 0;
                q.push({0, x});
                res[x] = {x, 0};
                ret[x] = 0;
            }

            while(!q.empty())
            {
				auto [d, a] = q.top();
                q.pop();

                if(dist[a] < d) continue;

 				for(auto [b, d1]: graph.GetAllEdges(a))
                {
                    // if(d1 + d < dist[b])
                    if(d1 + d < dist[b] || (d1 + d == dist[b] && res[b].first > res[a].first))
                    {
                        dist[b] = d1 + d;
                        q.push({d1 + d, b});

                        ret[b] = d1 + d;
                        res[b] = { res[a].first, d1 + d };
                    }
                }
            }
            return ret;
        }

        DistanceType Query(NodeType u, NodeType v)
        {
            DistanceType mindis=inf;
            if(L[u].size()>L[v].size()) swap(u,v);
            for(auto [x,w]: L[u])
                qtable[x]=w;
            for(auto [x,w]: L[v])
                if(qtable[x] != inf)
                    mindis=min(mindis,qtable[x]+w);
            for(auto [x,w]: L[u])
                qtable[x]=inf;
            return mindis;
        }

    
};