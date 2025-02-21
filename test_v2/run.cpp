#include "../HL/HL_basic.hpp"
#include "../ADO/ADO_basic.hpp"
#include "../ADO/ADO_fixed_Ak.hpp"
#include "../HL_ADO/mix_basic.hpp"
#include "../HL_ADO/mix_v2.hpp"
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
using std::map;
using std::string;
using std::cerr;
using Clock = std::chrono::high_resolution_clock;


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

map <string, int> graph_offset;
map <string, string> graph_type;
int seed;

void Load_Parameters()
{
    graph_offset["Toronto"] = 1;
    graph_offset["Movielens"] = 1;
    graph_offset["DBLP"] = 1;
    graph_offset["Dbpedia"] = 0;
    graph_offset["LinkedMDB"] = 0;
    graph_offset["LUBM-50K"] = 0;
    graph_offset["LUBM-500K"] = 0;
    graph_offset["LUBM-5M"] = 0;
    graph_offset["Mondial"] = 0;
    graph_offset["OpenCyc"] = 0;
    graph_offset["YAGO"] = 0;
    graph_offset["example"] = 0;

    graph_type["0"] = graph_type["UW"] = "UW";
    graph_type["1"] = graph_type["IW"] = "IW";

}

template<typename NodeType, typename DistanceType>
void RunHL_basic(AdjacencyList<NodeType, DistanceType> &adjList, FileManager &fop, string outfile, string resultname)
{
    int N = adjList.GetSize();

    cerr << "Start Build HL\n";

    auto start = Clock::now();
    HL_basic<NodeType, DistanceType> oracle(adjList);
    auto dur = Clock::now() - start;

    std::string create_time = "HL_basic generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "HL Size is " + std::to_string(oracle.LabelSize));
    
    cerr << "Start Output HL\n";
    oracle.Output_Label(outfile);

    cerr << "HL_basic End\n";
}

void Build_HL(string graphname, string type)
{
    // usage: run <graph-name> <graph-type>

    if(graph_offset.find(graphname) == graph_offset.end() || graph_type.find(type) == graph_type.end())
    {
        cerr << "usage: run <graph-name> <graph-type>\n";
        exit(0); 
    }

    //cerr << graph_type[type] << "\n";

    string graph_file = "../data/" + graphname;
    FileManager fop(graphname + "/" + graph_type[type]);
    
    string result_name = "HL_Info.txt";
    fop.clearFile(result_name);

    if(graph_type[type] == "UW")
    {
        LoadGraph <int, int> g(graph_file, 0, graph_offset[graphname]);
        RunHL_basic(g.adjList, fop, graph_file + "/HL_UW.txt", result_name);
    }
    else
    {
        LoadGraph <int, double> g(graph_file, 1, graph_offset[graphname]);
        RunHL_basic(g.adjList, fop, graph_file + "/HL_IW.txt", result_name);
    }
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> Get_Exact_Distance(AdjacencyList<NodeType, DistanceType> &adjList, string filepre, vector <pair <NodeType,NodeType> > &qry, FileManager &fop, string resultname)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);

    cerr << "Start Read HL\n";
    auto start = Clock::now();
    HL_basic<NodeType, DistanceType> oracle(adjList, filepre);
    auto dur = Clock::now() - start;

    cerr << "Start Query\n";
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average HL query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.\n";
    fop.writeToFile(resultname, query_time);
    cerr << "Get_Exact_Distance End\n";
    return ans;
}

template<typename NodeType, typename DistanceType>
std::pair <vector<DistanceType>,vector<DistanceType> > RunADO_v2(AdjacencyList<NodeType, DistanceType> &adjList, int K, vector <pair <NodeType,NodeType> > &qry, FileManager &fop, string resultname)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m),ans_(m);

    int N = adjList.GetSize();

    auto start = Clock::now();
    ADO_fixed_1<NodeType, DistanceType> oracle(adjList, K);
    oracle.setseed(seed);

    auto dur = Clock::now() - start;

    std::string create_time = "ADO_v2 generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "ADO Size is " + std::to_string(oracle.LabelSize));

    
    start = Clock::now();
    int gcnt = 0;
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query_Length(u, v);
        ans[i] = c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average ADO_Distance query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);


    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = ans[i];
        auto [c_, vex] = oracle.Query_Path(u, v);
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
        
        ans_[i] = c_;
    }

    dur = Clock::now() - start;

    std::cerr << "Great distance ratio is " << 1.0*gcnt/m << "\n";
    
    query_time = "Average ADO_Path query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);

    return std::make_pair(ans, ans_);
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunMix_Basic(AdjacencyList<NodeType, DistanceType> &adjList,vector <pair <NodeType,NodeType> > &qry, FileManager &fop, double ratio, string resultname)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    int N = adjList.GetSize();

    auto start = Clock::now();
    Mix_basic<NodeType, DistanceType> oracle(adjList, ceil(pow(N, ratio)) + 1);
    auto dur = Clock::now() - start;

    string Info = "Mix_Basic_n^{" + std::to_string(ratio) + "}";

    std::string create_time = Info + " generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, Info + " size is " + std::to_string(oracle.LabelSize));
    
    
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average " + Info + " query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);
    return ans;
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunMix_v2(AdjacencyList<NodeType, DistanceType> &adjList,int K,vector <pair <NodeType,NodeType> > &qry, FileManager &fop, string resultname)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    int N = adjList.GetSize();

    auto start = Clock::now();
    Mix_v2<NodeType, DistanceType> oracle(adjList, K);
    oracle.setseed(seed);
    auto dur = Clock::now() - start;

    std::string create_time = "Mix_v2 generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "Mix_v2 Size is " + std::to_string(oracle.LabelSize));
    
    
    start = Clock::now();
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }
    dur = Clock::now() - start;

    std::string query_time = "Average Mix_v2 query time: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()/m) + " ms.";
    fop.writeToFile(resultname, query_time);
    return ans;
}

template <typename DistanceType>
void Compare_Appro(vector <DistanceType> &anshl, vector <DistanceType> &ansap, FileManager &fop, int K, string Info, string resultname, int is_next_line)
{
    int qn=anshl.size();
    DistanceType inf = std::numeric_limits<DistanceType>::max();
    
    double ratio = 0;
    int cntc = 0;
    for(int i = 0; i < qn; i++)
    {
        DistanceType x = anshl[i], y = ansap[i];
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
            if(1.0*y/x > 2*K - 1)
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
    std::string info = "Average ratio is " + std::to_string(ratio) 
                        + " and correctness rate is " + std::to_string(1.0*cntc/qn);
    fop.writeToFile(resultname, Info + ": " + info); 
    if(is_next_line) fop.writeToFile(resultname, "");
}


void Run_Query(int argc,char* argv[])
{
    //usage: run <seed> <graph-name> <graph-type> <qry-number> <K> (<ratio>) : n^ratio
    if(argc != 6 && argc != 7)
    {
        cerr << "usage: run <seed> <graph-name> <graph-type> <qry-number> <K> (<ratio>)\n";
        exit(0); 
    }

    seed = atoi(argv[1]);
    string graphname = argv[2];
    string type = argv[3];
    const int QM = atoi(argv[4]);
    const int K = atoi(argv[5]);
    double ratio = -1;
    if(argc == 7) ratio = atof(argv[6]);

    if(graph_offset.find(graphname) == graph_offset.end() || graph_type.find(type) == graph_type.end())
    {
        cerr << "usage: run <seed> <graph-name> <graph-type> <qry-number> <K> (<ratio>)\n";
        exit(0); 
    }

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
    
    string graph_file = "../data/" + graphname;

    if(!fs::exists(graph_file + "/HL_" + graph_type[type]+ ".txt"))
        Build_HL(graphname, type);

    FileManager fop(graphname + "/" + graph_type[type]);
    
    string result_name = "Result_1.txt";
    fop.clearFile(result_name);

    auto read_hl_info = [&]()
    {
        string hl_info_file = graphname + "/" + graph_type[type] + "/HL_Info.txt";

        ifstream inputFile(hl_info_file);
        if (!inputFile.is_open()) 
            cerr << "Can't open file: " << hl_info_file << endl;
        inputFile.tie(nullptr);

        string Line;
        while(getline(inputFile, Line))
            fop.writeToFile(result_name, Line);
    };
    read_hl_info();

    auto GoTest = [&]<typename EdgeType>(const string& hl_suffix) 
    {
        LoadGraph<int, EdgeType> g(graph_file, graph_type[type] == "UW" ? 0 : 1, graph_offset[graphname]);
        N = g.N;
        auto qry = geneqry();
        auto ans = Get_Exact_Distance(g.adjList, graph_file + "/HL_" + hl_suffix + ".txt", qry, fop, result_name);

        auto [ado_d, ado_p] = RunADO_v2(g.adjList, K, qry, fop, result_name);
        Compare_Appro(ans, ado_d, fop, K, "ADO only distance", result_name, 0);
        Compare_Appro(ans, ado_p, fop, K, "ADO path", result_name, 1);

        auto mix_1 = RunMix_Basic(g.adjList, qry, fop, 0.5, result_name);
        Compare_Appro(ans, mix_1, fop, 2, "Mix sqrt", result_name, 1);

        auto mix_2 = RunMix_Basic(g.adjList, qry, fop, 0.4, result_name);
        Compare_Appro(ans, mix_2, fop, 2, "Mix ^0.4", result_name, 1);

        auto mix_v2 = RunMix_v2(g.adjList, K, qry, fop, result_name);
        Compare_Appro(ans, mix_v2, fop, K, "Mix_v2", result_name, 1);
    };

    if (graph_type[type] == "UW")
        GoTest.template operator()<int>("UW");
    else
        GoTest.template operator()<double>("IW");
    
}

int main(int argc, char* argv[])
{
    Load_Parameters();
    if(argc == 3)
        Build_HL(argv[1], argv[2]);
    else
        Run_Query(argc, argv);
    return 0;
}