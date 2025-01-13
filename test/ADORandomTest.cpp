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
vector<DistanceType> RunADO_basic(AdjacencyList<NodeType, DistanceType> adjList, int K, vector <pair <NodeType,NodeType> > qry)
{
    int N = adjList.GetSize();

    auto start = Clock::now();
    ADO_basic<NodeType, DistanceType> oracle(adjList, K);
    auto dur = Clock::now() - start;
    std::cout << "ADO_basic generated in: " << std::chrono::duration<double, std::milli>(dur).count() << " ms."<< std::endl;

    int m=qry.size();
    std::vector <DistanceType> ans(m);
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        // assert(c <= (2 * K - 1) * opt);
        ans[i]=c;
    }
    dur = Clock::now() - start;
    std::cout << "Average ADO query time: " <<  std::chrono::duration<double, std::milli>(dur).count()/m << " ms." << std::endl;
    return ans;
}


int main(int argc, char* argv[])
{
    if(argc <= 7)
    {
        std::cout << "usage: RandomGraphTest K <seed> N M W IN <qry-number>" << std::endl;
        std::cout << "For unweighted graph, use W = 0, otherwith the edge value will be randomly in [0,W]" << std::endl;
        std::cout << "For big deg & small diameter, use IN ~ sqrt(N), otherwise for similar roadgraph, use IN ~ O(1)" << std::endl;
        return -1;
    }

    const int K = atoi(argv[1]);
    const int seed = atoi(argv[2]);
    const int N = atoi(argv[3]);
    const int M = atoi(argv[4]);
    const double W = strtof(argv[5], nullptr);
    const int IN = atoi(argv[6]);
    const int QM = atoi(argv[7]);

    auto geneqry = [&]()
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> rnd(0, N-1);
        vector <pair <int,int> > qry(QM);
        for(int i=0;i<QM;i++)
            qry.push_back({rnd(gen),rnd(gen)});
        return qry;
    };

    auto start = Clock::now();

    if(W == 0)
    {
        auto adjList = GenGraph<int>(N, M, 1, 1.0*IN/N, seed);
        std::chrono::duration<double, std::milli> dur = Clock::now() - start;
        std::cout << "Graph of " << adjList.GetSize() << " vertices and " << adjList.GetEdgeCount() << " edges. Generated in: " << dur.count() << " ms." << std::endl;

        RunADO_basic(adjList, 3, geneqry());
    }
    else
    {
        auto adjList = GenGraph<double>(N, M , W, 1.0*IN/N, seed);
        std::chrono::duration<double, std::milli> dur = Clock::now() - start;
        std::cout << "Graph of " << adjList.GetSize() << " vertices and " << adjList.GetEdgeCount() << " edges. Generated in: " << dur.count() << " ms." << std::endl;

        RunADO_basic(adjList, 3, geneqry());
    }
    
    return 0;
}
