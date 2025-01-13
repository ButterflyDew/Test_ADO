#include "../HL/HL_basic.hpp"
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
#include <fstream>
#include <filesystem>
#include <map>
using std::pair; 
using std::vector;
using Clock = std::chrono::high_resolution_clock;
#define OnlyADO 1
#define OnlyHL 2
#define Nop 0

int ban;

namespace fs = std::filesystem;
class FileManager 
{
private:
    std::string filepath;

    void createDirectoryIfNotExist() 
    {
        if (!fs::exists(filepath)) 
        {
            std::cout << "Directory does not exist. Creating it..." << std::endl;
            fs::create_directories(filepath);
        }
    }

public:
    FileManager(const std::string& path) : filepath(path) 
    {
        createDirectoryIfNotExist();
    }

    void writeToFile(const std::string& filename, const std::string& info) 
    {
        std::string file = filepath + "/" + filename;

        std::ofstream outFile(file, std::ios::app);
        if (!outFile) 
        {
            std::cerr << "Failed to open the file." << std::endl;
            return;
        }

        outFile << info << std::endl;
    }

    void clearFile(const std::string& filename) {
        std::string file = filepath + "/" + filename;
        std::ofstream outFile(file, std::ios::trunc); 
    }
};


template<typename NodeType, typename DistanceType>
vector<DistanceType> RunHL_basic(AdjacencyList<NodeType, DistanceType> &adjList,vector <pair <NodeType,NodeType> > &qry, FileManager &fop)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    if(ban&1) 
    {
        std::cerr<<"Not Run HL"<<std::endl;
        return ans;
    }
    int N = adjList.GetSize();

    auto start = Clock::now();
    HL_basic<NodeType, DistanceType> oracle(adjList);
    auto dur = Clock::now() - start;

    std::string create_time = "HL_basic generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile("Result.txt", create_time);
    fop.writeToFile("Result.txt", "HL Size is " + std::to_string(oracle.LabelSize));
    
    
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average HL query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile("Result.txt", query_time);
    return ans;
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunADO_basic(AdjacencyList<NodeType, DistanceType> &adjList, int K, vector <pair <NodeType,NodeType> > &qry, FileManager &fop)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    if(ban>>1&1) 
    {
        std::cerr<<"Not Run ADO"<<std::endl;
        return ans;
    }
    int N = adjList.GetSize();

    auto start = Clock::now();
    ADO_basic<NodeType, DistanceType> oracle(adjList, K);
    auto dur = Clock::now() - start;

    std::string create_time = "ADO_basic generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile("Result.txt", create_time);
    fop.writeToFile("Result.txt", "ADO Size is " + std::to_string(oracle.LabelSize));

    
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

    std::string query_time = "Average ADO query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile("Result.txt", query_time);

    return ans;
}

template <typename DistanceType>
void Compare_Appro(vector <DistanceType> &anshl, vector <DistanceType> &ansado, FileManager &fop, int K)
{
    if(ban!=0) return;
    int qn=anshl.size();
    DistanceType inf = std::numeric_limits<DistanceType>::max();
    double ratio = 0;
    for(int i = 0; i < qn; i++)
    {
        DistanceType x = anshl[i], y = ansado[i];
        if(x == inf || y == inf)
        {
            if(x != inf || y != inf) 
            {
                std::cerr << x << " " << y << std::endl;
                std::cerr << "Error1" << std::endl;
                exit(0);
            }
            ratio += 1;
        }
        else if(x!=0)
        {
            if(1.0*y/x > 2*K - 1)
            {
                std::cerr << x << " " << y << std::endl;
                std::cerr << "Error2" << std::endl;
                exit(0);
            }
            ratio += 1.0*y/x;
        }    
    }
    ratio/=1.0*qn;
    std::string info = "Average ratio is " + std::to_string(ratio);
    fop.writeToFile("Result.txt", info); 
}
std::map <std::string, int> offs;
int main(int argc, char* argv[])
{
    if(argc != 6)
    {
        std::cout << "usage: Compare <seed> <graph-name> <graph-type> <qry-number> <K>" << std::endl;
        return -1;
    }

    const int seed = atoi(argv[1]);
    std::string graphname = argv[2];
    const int type = atoi(argv[3]);
    const int QM = atoi(argv[4]);
    const int K = atoi(argv[5]);
    offs["Toronto"] = offs["Movielens"] = offs["DBLP"] = 1;
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

    FileManager fop(graphname+(type == 0 ?"_UW":"_IW"));
    ban=OnlyADO;
    fop.clearFile("Result.txt");
    if(type==0)
    {
        LoadGraph <int,int> g(filename, 0, offs[graphname]);
        N = g.N;
        auto qry = geneqry();
        auto anshl = RunHL_basic(g.adjList, qry, fop);
        auto ansado = RunADO_basic(g.adjList, K, qry, fop);
        Compare_Appro(anshl, ansado, fop, K);
    }
    else
    {
        LoadGraph <int,double> g(filename, 1, offs[graphname]);
        N = g.N;
        auto qry = geneqry();
        auto anshl = RunHL_basic(g.adjList, qry, fop);
        auto ansado = RunADO_basic(g.adjList, K, qry, fop);
        Compare_Appro(anshl, ansado, fop, K);
    }
    
    return 0;
}
