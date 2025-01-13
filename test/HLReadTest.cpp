#include "../HL/HL_basic.hpp"
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
vector<DistanceType> RunHL_basic(AdjacencyList<NodeType, DistanceType> &adjList,vector <pair <NodeType,NodeType> > &qry)
{
    int N = adjList.GetSize();

    auto start = Clock::now();
    HL_basic<NodeType, DistanceType> oracle(adjList);
    auto dur = Clock::now() - start;
    std::cout << "HL_basic generated in: " << std::chrono::duration<double, std::milli>(dur).count() << " ms./"
    << std::chrono::duration<double>(dur).count() << "s." <<std::endl;

    int m=qry.size();
    std::vector <DistanceType> ans(m);
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;
    std::cout << "Average ADO query time: " <<  std::chrono::duration<double, std::milli>(dur).count()/m << " ms." << std::endl;
    return ans;
}


int main(int argc, char* argv[])
{
    if(argc <= 3)
    {
        std::cout << "usage: HL_Read_Test <seed> <graph-name> <graph-type> <qry-number>" << std::endl;
        return -1;
    }

    const int seed = atoi(argv[1]);
    std::string graphname = argv[2];
    const int type = atoi(argv[3]);
    const int QM = atoi(argv[4]);
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
        auto qry = geneqry();
        RunHL_basic(g.adjList, qry);
    }
    else
    {
        LoadGraph <int,double> g(filename, 1);
        N = g.N;
        auto qry = geneqry();
        RunHL_basic(g.adjList, qry);
    }
    
    return 0;
}
