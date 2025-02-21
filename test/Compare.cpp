#include "../HL/HL_basic.hpp"
#include "../ADO/ADO_basic.hpp"
#include "../ADO/ADO_fixed_Ak.hpp"
#include "../HL_ADO/mix.hpp"
#include "../HL_ADO/mix_fixed.hpp"
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
int seed;
std::string resultname;


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
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "HL Size is " + std::to_string(oracle.LabelSize));
    
    
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
    fop.writeToFile(resultname, query_time);
    return ans;
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunMix_Basic(AdjacencyList<NodeType, DistanceType> &adjList,vector <pair <NodeType,NodeType> > &qry, FileManager &fop)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    int N = adjList.GetSize();

    auto start = Clock::now();
    Mix_basic<NodeType, DistanceType> oracle(adjList, ceil(pow(N, 0.4)) + 1);
    auto dur = Clock::now() - start;

    std::string create_time = "Mix_Basic generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "Size is " + std::to_string(oracle.LabelSize));
    
    
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average Mix query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);
    return ans;
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunMix_fixed(AdjacencyList<NodeType, DistanceType> &adjList,int K,vector <pair <NodeType,NodeType> > &qry, FileManager &fop)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    int N = adjList.GetSize();

    auto start = Clock::now();
    Mix_fixed<NodeType, DistanceType> oracle(adjList, K);
    oracle.setseed(seed);
    auto dur = Clock::now() - start;

    std::string create_time = "Mix_fixed generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "Mix_fixed Size is " + std::to_string(oracle.LabelSize));
    
    
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average Mix query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);
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
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "ADO Size is " + std::to_string(oracle.LabelSize));

    
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        // assert(c <= (2 * K - 1) * opt);
        //std::cerr << "In basic, The Aprov. distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average ADO query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);

    return ans;
}

template<typename NodeType, typename DistanceType>
std::pair <vector<DistanceType>,vector<DistanceType> > RunADO_fixed_1(AdjacencyList<NodeType, DistanceType> &adjList, int K, vector <pair <NodeType,NodeType> > &qry, FileManager &fop)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m),ans_(m);
    if(ban>>1&1) 
    {
        std::cerr<<"Not Run ADO"<<std::endl;
        return std::make_pair(ans, ans_);
    }
    int N = adjList.GetSize();

    auto start = Clock::now();
    ADO_fixed_1<NodeType, DistanceType> oracle(adjList, K);
    oracle.setseed(seed);

    auto dur = Clock::now() - start;

    std::string create_time = "ADO_fixed_Ak generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "ADO Size is " + std::to_string(oracle.LabelSize));

    
    start = Clock::now();
    int gcnt = 0;
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query_Length(u, v);
        auto [c_, vex] = oracle.Query_Path(u, v);
        //assert(c <= (2 * K - 1) * opt);
        assert(c_ <= c);

        if(c_ != c) ++gcnt;

        if(!vex.empty())
        {
            int vlen = vex.size();
            DistanceType sumw = 0;
            for(int j = 1; j < vlen; j++) 
            {
                assert(adjList.GetEdgeuv(vex[j - 1], vex[j]) != std::numeric_limits<DistanceType>::max());
                sumw += adjList.GetEdgeuv(vex[j - 1], vex[j]);
            }
            if(abs(sumw-c_) > 1e-9)
                std::cerr << std::fixed << std::setprecision(10) << sumw << " " << c_ << "\n";
            assert(abs(sumw-c_) <= 1e-9);
        }
        
        // std::cerr << "In f1, The Aprov. distance of " << u << " and " << v << " is " << c << "\n";
        // std::cerr << "c_ is " << c_ << "\n" << "Path is:";
        // for(auto x: vex) std::cerr << x << " ";
        // std::cerr << "\n\n";
        ans[i] = c;
        ans_[i] = c_;
    }
    std::cerr << "Great distance ratio is " << 1.0*gcnt/m << "\n";
    dur = Clock::now() - start;

    std::string query_time = "Average ADO query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);

    return std::make_pair(ans, ans_);
}

template <typename DistanceType>
void Compare_Appro(vector <DistanceType> &anshl, vector <vector <DistanceType> > &ansado, FileManager &fop, vector <int> K)
{
    if(ban!=0) return;
    int qn=anshl.size();
    DistanceType inf = std::numeric_limits<DistanceType>::max();
    int cur = 0;
    for(auto ado: ansado)
    {
        int k_now = K[cur];
        ++cur;
        double ratio = 0;
        int cntc = 0;
        for(int i = 0; i < qn; i++)
        {
            DistanceType x = anshl[i], y = ado[i];
            if(abs(x-y) < 1e-9) ++cntc; 
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
                if(1.0*y/x > 2*k_now - 1)
                {
                    std::cerr << x << " " << y << std::endl;
                    std::cerr << "Error2" << std::endl;
                    exit(0);
                }
                ratio += 1.0*y/x;
            }    
            else
                ratio += 1;
        }
        ratio/=1.0*qn;
        std::string info = "Average ratio " + std::to_string(cur) + " is " + std::to_string(ratio) 
                         + " and correctness rate is " + std::to_string(1.0*cntc/qn);
        fop.writeToFile(resultname, info); 
    }
}
std::map <std::string, int> offs;
int main(int argc, char* argv[])
{
    if(argc != 6)
    {
        std::cout << "usage: Compare <seed> <graph-name> <graph-type> <qry-number> <K>" << std::endl;
        return -1;
    }

    seed = atoi(argv[1]);
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
    FileManager fop(graphname+(type == 0 ?"_UW":"_IW") + "_K=" + std::to_string(K));
    
    
    ban = Nop;
    resultname = "Result_fixed5.txt";
    
    fop.clearFile(resultname);
    if(type==0)
    {
        LoadGraph <int,int> g(filename, 0, offs[graphname]);
        N = g.N;
        auto qry = geneqry();
        auto anshl = RunHL_basic(g.adjList, qry, fop);
        auto ansado = RunADO_basic(g.adjList, K, qry, fop);
        auto ansado_f1 = RunADO_fixed_1(g.adjList, K, qry, fop);
        auto ansmix = RunMix_Basic(g.adjList, qry, fop);
        //auto ansmix = RunMix_fixed(g.adjList, K, qry, fop);
        vector <vector <int> > ado;
        vector <int> Ki;
        ado.push_back(ansado), ado.push_back(ansado_f1.first), ado.push_back(ansado_f1.second), ado.push_back(ansmix);
        Ki.push_back(K), Ki.push_back(K), Ki.push_back(K), Ki.push_back(K);
        Compare_Appro(anshl, ado, fop, Ki);
    }
    else
    {
        LoadGraph <int,double> g(filename, 1, offs[graphname]);
        N = g.N;
        auto qry = geneqry();
        auto anshl = RunHL_basic(g.adjList, qry, fop);
        auto ansado = RunADO_basic(g.adjList, K, qry, fop);
        auto ansado_f1 = RunADO_fixed_1(g.adjList, K, qry, fop);
        auto ansmix = RunMix_Basic(g.adjList, qry, fop);
        //auto ansmix = RunMix_fixed(g.adjList, K, qry, fop);
        vector <vector <double> > ado;
        vector <int> Ki;
        ado.push_back(ansado), ado.push_back(ansado_f1.first), ado.push_back(ansado_f1.second), ado.push_back(ansmix);
        Ki.push_back(K), Ki.push_back(K), Ki.push_back(K), Ki.push_back(K);
        Compare_Appro(anshl, ado, fop, Ki);
    }
    
    return 0;
}
