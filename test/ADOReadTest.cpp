#include "../ADO/ADO_basic.hpp"
#include "../AdjacencyList.hpp"
#include "../Utilities.hpp"
#include "Load_Graph.hpp"

#include <cassert>
#include <iostream>
#include <limits>
#include <random>
#include <chrono>
#include <algorithm>
#include <variant>
using std::pair; 
using std::vector;
using Clock = std::chrono::high_resolution_clock;

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunADO_basic(AdjacencyList<NodeType, DistanceType> &adjList, int K, vector <pair <NodeType,NodeType> > qry)
{
    int N = adjList.GetSize();

    auto start = Clock::now();
    ADO_basic<NodeType, DistanceType> oracle(adjList, K);
    auto dur = Clock::now() - start;
    std::cout << "ADO_basic generated in: " << std::chrono::duration<double, std::milli>(dur).count() << " ms./"
    << std::chrono::duration<double>(dur).count() << "s." <<std::endl;
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        // assert(c <= (2 * K - 1) * opt);
        //std::cerr << "The Aprov. distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;
    std::cout << "Average ADO query time: " <<  std::chrono::duration<double, std::milli>(dur).count()/m << " ms." << std::endl;
    return ans;
}


int main(int argc, char* argv[])
{
    if(argc <= 4)
    {
        std::cout << "usage: ADO_Read_Test K <seed> <graph-name> <graph-type> <qry-number>" << std::endl;
        return -1;
    }

    const int K = atoi(argv[1]);
    const int seed = atoi(argv[2]);
    std::string graphname = argv[3];
    const int type = atoi(argv[4]);
    const int QM = atoi(argv[5]);
    int N;
    auto geneqry = [&]()
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> rnd(0, N-1);
        vector <pair <int,int> > qry;
        for(int i=0;i<QM;i++)
            qry.push_back({rnd(gen),rnd(gen)});
        return qry;
    };

    std::string filename = "../data/"+graphname;
    if(type==0)
    {
        LoadGraph <int,int> g(filename, 0);
        N = g.N;
        RunADO_basic(g.adjList, 3, geneqry());
    }
    else
    {
        LoadGraph <int,double> g(filename, 1);
        N = g.N;
        RunADO_basic(g.adjList, 3, geneqry());
    }
    
    return 0;
}
