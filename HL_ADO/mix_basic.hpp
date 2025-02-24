#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <limits>
#include "../AdjacencyList.hpp"

using std::vector;
using std::pair;
using std::sort;
using std::min;
using std::swap;

template<typename NodeType, typename DistanceType>
class Mix_basic
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    private:
        AdjList &graph;
        vector<vector<pair<NodeType, DistanceType>>> L;
        vector<DistanceType> qtable;
        DistanceType inf;

    public:
        long long LabelSize;
        Mix_basic(AdjList &lis, int Bl): graph(lis)
        {
            LabelSize = 0;
            std::cerr << "Start Mix!" << std::endl;
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

            vector <NodeType> isr(n);
            vector <DistanceType> d(n,inf);

            vector <int> Ai;

            for(int i=0;i<Bl;i++)
            {
                vector <NodeType> inq;
                NodeType rt=vid[i];
                Ai.push_back(vid[i]);
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
                            if(d[u]+w<d[v])
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

            int psize = 0;
            for(int i = 0; i < n; i++) psize += L[i].size();
            
            std::cerr << "The first Bl vextices Label Size: " << psize << "\n";

            auto pi = SSSPTree(Ai);
            for(int i=Bl;i<n;i++)
            {

                vector <NodeType> inq;
                NodeType rt=vid[i];
                Ai.push_back(vid[i]);
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

            std::cerr << "Other vextices Label Size: " << LabelSize - psize << "\n";
            std::cerr << "Mix Contruct Done!" << std::endl;
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

        bool Query_Greater(NodeType u, DistanceType d)
        {
            for(auto [v, w]: L[u]) if(qtable[v] <= d - w) return false;
            return true;
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