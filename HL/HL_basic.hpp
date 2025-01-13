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
class HL_basic
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    public:
        long long LabelSize;
        HL_basic(AdjList &lis): graph(lis)
        {
            LabelSize = 0;
            std::cerr << "Start Construct HL!" << std::endl;
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

            
            vector <DistanceType> d(n,inf);

            int up = (n+9)/10, cur = 1;

            for(int i=0;i<n;i++)
            {
                if(i == up*cur)
                {
                    std::cerr << cur << "0% has done!" << std::endl;
                    ++cur;
                }

                vector <NodeType> inq;
                NodeType rt=vid[i];
                d[rt]=0;
                std::priority_queue <pair <DistanceType, NodeType>, vector <pair <DistanceType, NodeType> >, std::greater <pair <DistanceType, NodeType>> > q;
                q.push({d[rt],rt});
                while(!q.empty())
                {
                    auto [dis, u]=q.top();
                    q.pop();
                    if(d[u]<dis) continue;
                    inq.push_back(u);
                    if(d[u] < Query(rt, u))
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
            }

            for(int i = 0; i < n; i++) LabelSize += L[i].size();
            std::cerr << "HL Contruct Done!" << std::endl;
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

    private:
        AdjList &graph;
        vector<vector<pair<NodeType, DistanceType>>> L;
        vector<DistanceType> qtable;
        DistanceType inf;
};