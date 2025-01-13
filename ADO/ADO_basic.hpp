#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>
#include <unordered_map>

#include "../AdjacencyList.hpp"

template<typename NodeType, typename DistanceType>
class ADO_basic
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    public:
        long long LabelSize;
        ADO_basic(AdjList &lst, int K) : graph(lst), K(K)
        {
            LabelSize = 0;

            dmx = std::numeric_limits<DistanceType>::max(); 
            std::cerr << "Start Construct ADO!" << std::endl;
            A.assign(K + 1, std::vector<NodeType>());

            // used for fast lookup for calculating A[i] - A[i + 1]
            AC.assign(K + 1, std::vector<bool>(graph.GetSize(), false));

            B.assign(graph.GetSize(), std::unordered_map<NodeType, DistanceType>());

            P.assign(K + 1, std::vector<std::pair<NodeType, DistanceType>>(graph.GetSize(), std::make_pair(-1, dmx)));

            f.assign(graph.GetSize(), 0);

            dtmp.assign(graph.GetSize(), dmx);

            for(NodeType i = 0; i < graph.GetSize(); i++)
            {
                f[i]=i;
                A[0].push_back(i);
                AC[0][i] = true;
            }

            for(NodeType u = 0; u < graph.GetSize(); u++)
                for(auto v: graph.GetAllEdges(u))
                    Merge(u,v.first);

            double prob = std::pow((double)graph.GetSize(), -1.00 / K);
            std::bernoulli_distribution dist(prob);

            std::random_device rd;
            std::mt19937 gen(rd()); // seed

            for(int i = 1; i < K; i++)
            {
                do {
                    for(NodeType x: A[i - 1])
                    {
                        if(dist(gen))
                        {
                            A[i].push_back(x);
                            AC[i][x] = true;
                        }
                    }
                } while(A[i].size() == 0); // We dont want to get an empty list before A[K] but can this end up in infinite loop ?
            }
            // The expected size of A_i -> n^{1-i/k}, 
            std::cerr << "Sample A_0 ~ A_k Done!" << std::endl;

            std::vector<std::unordered_map<NodeType, DistanceType>> C(graph.GetSize());
            for(int i = K - 1; i >= 0; i--)
            {
                //Run Full SSSP with Source A[i]
                std::vector<std::pair<NodeType, DistanceType>> T = SSSPTree(i);

                for(NodeType x = 0; x < graph.GetSize(); x++)
                {
                    // Find nearest to x in A[i]
  					auto [piv, div] = T[x];
                    if(div == P[i + 1][x].second) piv = P[i + 1][x].first;

                    P[i][x] = {piv, div};
                }

                for(NodeType x: A[i])
                {
                    // skip if A[i + 1] contains x
                    if(AC[i + 1][x]) continue;

                    // grow shortest path tree from x
                    C[x] = GenCluster(x, i);
                }

                std::cerr << "The info of A[" << i << "] has Done!" << std::endl;
            }

            for(NodeType w = 0; w < graph.GetSize(); w++)
            {
                for(auto [v, d]: C[w])
                    B[v][w] = d;
            }

            for(NodeType u = 0; u < graph.GetSize(); u++)
                LabelSize += B[u].size();
            std::cerr << "ADO Contruct Done!" << std::endl;
        }

        // HACK: Do we need to care about P_i(v) = P_{i + 1}(v) in case of multiple witnesses, shouldn't priority_queue maintain that?
        std::vector<std::pair<NodeType, DistanceType>> SSSPTree(int i)
        {
            std::vector<DistanceType> dist(graph.GetSize(), dmx);

            std::priority_queue<std::pair<DistanceType, NodeType>, std::vector<std::pair<DistanceType, NodeType>>, std::greater<std::pair<DistanceType, NodeType>>> q;

            std::vector<std::pair<NodeType, DistanceType>> res(graph.GetSize(), {-1, dmx});

            for(NodeType x: A[i])
            {
                dist[x] = 0;
                q.push({0, x});
                res[x] = {x, 0};
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

                        res[b] = { res[a].first, d1 + d };
                    }
                }
            }
            return res;
        }

        std::unordered_map<NodeType, DistanceType> GenCluster(NodeType w, int pi)
        {
            // std::vector<DistanceType> dist(graph.GetSize(), std::numeric_limits<DistanceType>::max());
            // dist[w] = 0;
            dtmp[w] = 0;
            std::vector <NodeType> inq;
            std::priority_queue<std::pair<DistanceType, NodeType>, std::vector<std::pair<DistanceType, NodeType>>, std::greater<std::pair<DistanceType, NodeType>>> q;

            q.push({0, w});

            std::unordered_map<NodeType, DistanceType> res;
            while(!q.empty())
            {
				auto [d, a] = q.top();
                q.pop();

                if(dtmp[a] < d) continue;

                inq.push_back(a);
                res[a] = d;
 				for(auto [b, d1]: graph.GetAllEdges(a))
                {
                    if(d1 + d < dtmp[b] && d1 + d < P[pi + 1][b].second)
                    {
                        dtmp[b] = d1 + d;
                        q.push({d1 + d, b});
                    }
                }
            }
            for(auto a: inq) 
                dtmp[a] = dmx;
            return res;
        }

        DistanceType Query(NodeType u, NodeType v)
        {
            if(Find(u)!=Find(v)) return std::numeric_limits<DistanceType>::max();
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

    private:
        int K;
        AdjList& graph;
        std::vector<std::vector<NodeType>> A;
        std::vector<std::vector<bool>> AC;
        std::vector<std::vector<std::pair<NodeType, DistanceType>>> P;
        std::vector<std::unordered_map<NodeType, DistanceType>> B;
        std::vector<DistanceType> dtmp;
        std::vector<NodeType> f;
        DistanceType dmx;
        NodeType Find(NodeType x){return f[x]=f[x]==x?x:Find(f[x]);}
        void Merge(NodeType u,NodeType v)
        {
            if(Find(u)!=Find(v)) 
                f[Find(v)]=Find(u);
        }
};
