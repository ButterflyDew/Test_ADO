#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unordered_set>
#include <unordered_map>

#include "../AdjacencyList.hpp"

template<typename NodeType, typename DistanceType>
class ADO_fixed_1
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    private:
        int K;
        int seed;
        struct Cluster
        {
            NodeType x;
            DistanceType d;
            int par, dep;
            void Print(int cur){std::cerr << "In "<< cur << " x:" << x << " dis: "<< d << " par: " << par << " dep: " << dep << "\n";}
            Cluster(){}
            Cluster(NodeType x_): x(x_){d = 0, par = dep = -1;}
            Cluster(NodeType x_, DistanceType d_, int par_, int dep_): x(x_), d(d_), par(par_), dep(dep_) {}
            bool friend operator < (Cluster A, Cluster B)
            {
                return A.x < B.x;
            }
        };
        AdjList& graph;
        std::vector<std::vector<NodeType>> A;
        std::vector<std::vector<bool>> AC;
        std::vector<std::vector<std::pair<NodeType, DistanceType>>> P;
        std::vector <std::vector <std::pair<NodeType, DistanceType>>> B;
        std::vector <std::vector <Cluster> > C; 
        std::vector<DistanceType> dtmp;
        std::vector<int> deptmp;
        std::vector<NodeType> partmp;
        std::vector<NodeType> f;
        DistanceType dmx;
        NodeType Find(NodeType x){return f[x]=f[x]==x?x:Find(f[x]);}
        void Merge(NodeType u,NodeType v)
        {
            if(Find(u)!=Find(v)) 
                f[Find(v)]=Find(u);
        }


    public:
        long long LabelSize;

        void setseed(int Seed){seed=Seed;}

        ADO_fixed_1(AdjList &lst, int K) : graph(lst), K(K)
        {
            LabelSize = 0;

            dmx = std::numeric_limits<DistanceType>::max(); 
            std::cerr << "Start Construct ADO!" << std::endl;
            A.assign(K + 1, std::vector<NodeType>());

            // used for fast lookup for calculating A[i] - A[i + 1]
            AC.assign(K + 1, std::vector<bool>(graph.GetSize(), false));

            B.assign(graph.GetSize(), std::vector <std::pair <NodeType, DistanceType> >() );

            C.assign(graph.GetSize(), std::vector <Cluster> () );

            P.assign(K + 1, std::vector<std::pair<NodeType, DistanceType>>(graph.GetSize(), std::make_pair(-1, dmx)));

            f.assign(graph.GetSize(), 0);

            dtmp.assign(graph.GetSize(), dmx);
            deptmp.assign(graph.GetSize(), -1);
            partmp.assign(graph.GetSize(), -1);

            for(NodeType i = 0; i < graph.GetSize(); i++)
            {
                f[i]=i;
                A[0].push_back(i);
                AC[0][i] = true;
            }

            std::vector <NodeType> deg(graph.GetSize()), perm(graph.GetSize());

            for(NodeType u = 0; u < graph.GetSize(); u++)
            {
                perm[u]=u;
                for(auto v: graph.GetAllEdges(u))
                {
                    Merge(u,v.first);
                    ++deg[u], ++deg[v.first];
                }    
            }
            std::sort(perm.begin(), perm.end(), [&](auto x,auto y){return deg[x]>deg[y];});

            int Size = ceil(pow(graph.GetSize(), 1.0/K));
            perm.resize(Size);

            std::mt19937 gen(seed);
            shuffle(A[0].begin(), A[0].end(), gen);

            for(int i = 1; i < K; i++)
            {
                Size = ceil(pow(graph.GetSize(), 1.0 - 1.0*i/K));
                for(NodeType x: perm) A[i].push_back(x);
                int pos = A[i-1].size() - 1;
                while(A[i].size() < Size) A[i].push_back(A[i-1][pos--]);

                for(NodeType x: A[i]) AC[i][x] = true;
            }

            std::cerr << "Sample A_0 ~ A_k Done!" << std::endl;

            //std::vector<std::unordered_map<NodeType, DistanceType>> C(graph.GetSize());
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
                    GenCluster(x, i, C[x]);
                }

                std::cerr << "The info of A[" << i << "] has Done!" << std::endl;
            }

            for(NodeType w = 0; w < graph.GetSize(); w++)
            {
                //std::cerr << "w = " << w << "\n";
                //int curid = 0;
                for(auto Cx: C[w])
                {
                    // Cx.Print(curid);
                    // ++curid;
                    B[Cx.x].push_back({w, Cx.d});
                }
                // std::cerr << "\n";
            }

            for(NodeType u = 0; u < graph.GetSize(); u++)
            {
                sort(B[u].begin(), B[u].end());
                LabelSize += B[u].size();
            }    
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

        void GenCluster(NodeType w, int pi, std::vector <Cluster> &Cx)
        {
            // std::vector<DistanceType> dist(graph.GetSize(), std::numeric_limits<DistanceType>::max());
            // dist[w] = 0;
            dtmp[w] = 0;
            deptmp[w] = 0;
            partmp[w] = -1;
            std::vector <NodeType> inq;
            std::priority_queue<std::pair<DistanceType, NodeType>, std::vector<std::pair<DistanceType, NodeType>>, std::greater<std::pair<DistanceType, NodeType>>> q;

            q.push({0, w});
            while(!q.empty())
            {
				auto [d, a] = q.top();
                q.pop();

                if(dtmp[a] < d) continue;
                if(a != w)
                    deptmp[a] = deptmp[partmp[a]] + 1;
                Cx.emplace_back(a, dtmp[a], partmp[a], deptmp[a]);
                inq.push_back(a);
 				for(auto [b, d1]: graph.GetAllEdges(a))
                {
                    if(d1 + d < dtmp[b] && d1 + d < P[pi + 1][b].second)
                    {
                        dtmp[b] = d1 + d;
                        partmp[b] = a;
                        q.push({dtmp[b], b});
                    }
                }
            }

            sort(Cx.begin(), Cx.end());
            for(auto &x: Cx)
            {
                if(x.par == -1) continue;
                x.par = std::lower_bound(Cx.begin(), Cx.end(), Cluster(x.par)) - Cx.begin();
            }

            for(auto a: inq) 
                dtmp[a] = dmx;
        }

        NodeType Query_Mid(NodeType u, NodeType v)
        {
            auto is_in = [&](int v, int w)
            {
                auto it = std::lower_bound(B[v].begin(), B[v].end(), std::make_pair((NodeType)w, (DistanceType)-1));
                return it->first == w;
            };
            int i = 0;
            NodeType w = u;
            while(!is_in(v, w))
            {
                i++;
                std::swap(u, v);
                w = P[i][u].first;
            }
            return w;
        }

        DistanceType Query_Length(NodeType u, NodeType v)
        {
            if(Find(u)!=Find(v)) return std::numeric_limits<DistanceType>::max();
            
            NodeType w = Query_Mid(u, v);
            
            auto get_d = [&](int x, int w)
            {
                auto it = std::lower_bound(B[x].begin(), B[x].end(), std::make_pair((NodeType)w, (DistanceType)-1));
                return it->second;
            };

            return get_d(u, w) + get_d(v, w);
        }

        std::pair <DistanceType, std::vector <NodeType> > Query_Path(NodeType u, NodeType v)
        {
            if(Find(u)!=Find(v)) return {std::numeric_limits<DistanceType>::max(), std::vector <NodeType> ()};
            
            NodeType w = Query_Mid(u, v);
            
            int pu = std::lower_bound(C[w].begin(), C[w].end(), Cluster(u)) - C[w].begin();
            int pv = std::lower_bound(C[w].begin(), C[w].end(), Cluster(v)) - C[w].begin();

            //std::cerr << "w: " << w << " pu: " << pu << " pv: " << pv << "\n";

            if(C[w][pu].dep < C[w][pv].dep) std::swap(pu, pv);
            std::vector <NodeType> p1, p2;
            DistanceType sumw = C[w][pu].d + C[w][pv].d;
            while(C[w][pu].dep > C[w][pv].dep) p1.push_back(C[w][pu].x), pu = C[w][pu].par;
            while(pu != pv) 
            {
                p1.push_back(C[w][pu].x), pu = C[w][pu].par;
                p2.push_back(C[w][pv].x), pv = C[w][pv].par;
            }
            p1.push_back(C[w][pu].x);
            sumw -= C[w][pu].d*2;
            p1.insert(p1.end(), p2.rbegin(), p2.rend());
            return {sumw, p1};
        }
};
