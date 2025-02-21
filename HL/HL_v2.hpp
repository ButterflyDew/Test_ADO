#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <limits>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "../AdjacencyList.hpp"

using std::vector;
using std::pair;
using std::sort;
using std::min;
using std::swap;
using std::cerr;
using std::string;
using std::endl;
using std::istringstream;
using std::ifstream;
using std::ofstream;

template<typename NodeType, typename DistanceType>
class HL_v2
{
    using AdjList = AdjacencyList<NodeType, DistanceType>;

    public:
        long long LabelSize;
        HL_v2(AdjList &lis, string filepre): graph(lis)
        { 
            Read_Label(filepre);
        }
        HL_v2(AdjList &lis, int Sta): graph(lis), Status(Sta)
        {
            LabelSize = 0;
            std::cerr << "Start Construct HL!" << std::endl;
            n=graph.GetSize();
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

            vector <NodeType> isr(n), tag(n + 1), par(n + 1);
            
            vector <DistanceType> d(n,inf);

            int cur = 1;

            for(int i=0;i<n;i++)
            {
                if(i == (int)(1.0*n/10.0*cur))
                {
                    std::cerr << cur << "0% has done!" << std::endl;
                    ++cur;
                }

                vector <NodeType> inq;
                NodeType rt=vid[i];
                d[rt]=0;
                par[rt] = n;

                for(auto [v, w]: L[rt])
                    qtable[v] = w;

                std::priority_queue <pair <DistanceType, NodeType>, vector <pair <DistanceType, NodeType> >, std::greater <pair <DistanceType, NodeType>> > q;
                q.push({d[rt],rt});
                while(!q.empty())
                {
                    auto [dis, u]=q.top();
                    q.pop();
                    
                    if(d[u]<dis) continue;

                    inq.push_back(u);

                    if((Status&1)&&isr[u]) continue;

                    tag[u] = tag[par[u]]|isr[u];

                    if(((Status>>1&1) && tag[u] == 0) || ( (Status>>2&1) ? Query_Greater(u, d[u]) : d[u] < Query(rt, u) ) )
                    {
                        L[u].push_back({rt,d[u]});
                        for(auto [v,w]: graph.GetAllEdges(u))
                            if(d[u]+w<d[v])
                            {
                                d[v] = d[u] + w;
                                par[v] = u;
                                q.push({d[v],v});
                            }
                    }
                }
                for(auto u: inq) d[u] = inf, tag[u] = 0;
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

        void Output_Label(string filepre)
        {
            ofstream outputFile(filepre);
            outputFile << std::fixed << std::setprecision(15);
            outputFile << n << '\n';
            long long sum = 0;
            int per = 1;
            for(int i = 0; i < n; i++)
            {
                if(i == (int)(1.0*n/10.0*per))
                {
                    cerr << "Finish print "<< per*10 <<"%\n";
                    ++per;
                }
                for(auto [v,w]: L[i])
                    outputFile << v << " " << w <<  ",";
                outputFile << '\n';
                sum += L[i].size();
            }

            cerr << "The row of HBLL is " << n << '\n';
            cerr << "The sum of label is: " << sum << '\n';
            
            outputFile.close();
        }

        void Read_Label(string filepre)
        {
            clear_all();
            //cerr << filepre << "\n";
            ifstream inputFile(filepre);

            if (!inputFile.is_open()) 
                cerr << "Can't open file: " <<  filepre << endl;

            inputFile.tie(nullptr);

            inputFile >> n;
            L.assign(n, vector<pair <NodeType, DistanceType> >());
            qtable.assign(n, inf);
                
            int i = 0, per = 0;
            string Line;
            getline(inputFile, Line);
            while(getline(inputFile, Line))
            {
                if(i == (int)(1.0*n/10.0*per))
                {
                    cerr << "Finish read "<< per*10 <<"%\n";
                    ++per;
                }
                istringstream line(Line);
                string value;
                while(getline(line, value, ','))
                {
                    //cerr << value << "\n";
                    istringstream iss(value);
                    NodeType v; DistanceType w;
                    iss >> v >> w;
                    L[i].push_back({v,w});
                }
                ++i;
            }
            //cerr << "End and i is " << i << std::endl;
            inputFile.close();
        }

    private:
        AdjList &graph;
        vector<vector<pair<NodeType, DistanceType>>> L;
        vector<DistanceType> qtable;
        DistanceType inf;
        NodeType n;
        int Status;

        void clear_all()
        {
            n = 0;
            inf = std::numeric_limits<DistanceType>::max();
            L.clear();
            qtable.clear();
        }
};