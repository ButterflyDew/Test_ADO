#include "../HL/HL_basic.hpp"
#include "../HL/HL_v2.hpp"
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
void CK_HL_dijk(AdjacencyList<NodeType, DistanceType> &adjList, FileManager &fop, vector <NodeType> &qry, string outfile, string resultname)
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
    
    int m = qry.size();
    for(int i = 0; i < m; i++)
    {
        auto u = qry[i];
        auto dist = adjList.GetNearest(u);
        for(int j = 0; j < N; j++)
        {
            DistanceType c = oracle.Query(u, j);
            if(abs(c - dist[j]) > 1e-9)
            {
                cerr << "The exact distance of " << u << " and " << j << " are " << c << " and " << dist[j] << "\n";
                exit(0);
            }
            //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        }    
        
    }

    // cerr << "Start Output HL\n";
    // oracle.Output_Label(outfile);

    cerr << "Ck HL used dijk ok!\n";
}

template<typename NodeType, typename DistanceType>
vector<DistanceType> RunHL_basic(AdjacencyList<NodeType, DistanceType> &adjList, FileManager &fop, vector <pair <NodeType,NodeType> > &qry, string outfile, string resultname)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);

    int N = adjList.GetSize();

    cerr << "Start Build HL\n";

    auto start = Clock::now();
    HL_basic<NodeType, DistanceType> oracle(adjList);
    auto dur = Clock::now() - start;

    std::string create_time = "HL_basic generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "HL Size is " + std::to_string(oracle.LabelSize));
    
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }

    // cerr << "Start Output HL\n";
    // oracle.Output_Label(outfile);

    cerr << "HL_basic End\n";

    return ans;
}

template<typename NodeType, typename DistanceType>
std::vector <DistanceType> RunHL_v2(AdjacencyList<NodeType, DistanceType> &adjList, FileManager &fop, vector <pair <NodeType,NodeType> > &qry, string outfile, string resultname, string Info, int type)
{
    int m=qry.size();
    std::vector <DistanceType> ans(m);
    int N = adjList.GetSize();

    cerr << "Start Build HL_" + Info + "\n";

    auto start = Clock::now();
    HL_v2<NodeType, DistanceType> oracle(adjList, type);
    auto dur = Clock::now() - start;

    std::string create_time = "HL_" + Info + " generated in: " + std::to_string(std::chrono::duration<double, std::milli>(dur).count()) + "ms./"
    + std::to_string(std::chrono::duration<double>(dur).count()) + "s.";
    fop.writeToFile(resultname, create_time);
    fop.writeToFile(resultname, "HL_" + Info  + " Size is " + std::to_string(oracle.LabelSize));
    
    for(int i = 0; i < m; i++)
    {
        auto [u,v] = qry[i];
        DistanceType c = oracle.Query(u, v);
        //std::cerr << "The exact distance of " << u << " and " << v << " is " << c << "\n";
        ans[i]=c;
    }

    // cerr << "Start Output HL\n";
    // oracle.Output_Label(outfile);

    cerr << "HL_" + Info + " End\n";
    
    return ans;
}


template<typename NodeType, typename DistanceType>
void check(vector <DistanceType> &anshl, vector <DistanceType> &ansap, vector <pair <NodeType,NodeType> > &qry, string Info)
{
    cerr << "In " << Info << "\n";
    int qn=anshl.size();
    DistanceType inf = std::numeric_limits<DistanceType>::max();
    
    double ratio = 0;
    int cntc = 0;
    for(int i = 0; i < qn; i++)
    {
        DistanceType x = anshl[i], y = ansap[i];
        if(abs(x-y) > 1e-9)
        {
            auto [u, v] = qry[i];
            cerr << "The distance of " << u << " and " << v << " is: " << std::fixed << std::setprecision(15) << x << " and " << y << "\n";
            exit(0);
        }
    }
    cerr << Info << " ok\n";
}

void Build_HL(string graphname, string type, int QM)
{
    // usage: test_hl <graph-name> <graph-type> <qry-number>

    if(graph_offset.find(graphname) == graph_offset.end() || graph_type.find(type) == graph_type.end())
    {
        cerr << "usage: test_hl <graph-name> <graph-type> <qry-number>\n";
        exit(0); 
    }
    int N;

    auto geneNode = [&]()
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> rnd(0, N-1);
        vector <int> qry;
        for(int i=0;i<QM;i++)
            qry.push_back(rnd(gen));
        return qry;
    };

    
    auto geneqry = [&]()
    {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<> rnd(0, N-1);
        vector <pair <int,int> > qry;
        for(int i=0;i<QM;i++)
            qry.push_back({rnd(gen),rnd(gen)});
        return qry;
    };

    //cerr << graph_type[type] << "\n";

    string graph_file = "../data/" + graphname;
    FileManager fop(graphname + "/" + graph_type[type]);
    
    string result_name = "HL_Test_Info_dij_ckhl.txt";
    fop.clearFile(result_name);

    auto GoTest = [&]<typename EdgeType>(const string& hl_suffix) 
    {
        LoadGraph<int, EdgeType> g(graph_file, graph_type[type] == "UW" ? 0 : 1, graph_offset[graphname]);
        N = g.N;
        auto qry = geneqry();
        auto nlis = geneNode();
        string file = graph_file + "/HL_" + hl_suffix;
        CK_HL_dijk(g.adjList, fop, nlis, file + ".txt", result_name);
        // auto ans = RunHL_basic(g.adjList, fop, qry, file + ".txt", result_name);
        // auto a1 = RunHL_v2(g.adjList, fop, qry, file + "_1.txt", result_name, "1", 1);
        // auto a2 = RunHL_v2(g.adjList, fop, qry, file + "_2.txt", result_name, "2", 2);
        // auto a3 = RunHL_v2(g.adjList, fop, qry, file + "_3.txt", result_name, "3", 4);
        // auto aa = RunHL_v2(g.adjList, fop, qry, file + "_all.txt", result_name, "all", 7);
        // check(ans, a1, qry, "a1");
        // check(ans, a2, qry, "a2");
        // check(ans, a3, qry, "a3");
        // check(ans, aa, qry, "aa");
    };
    

    if (graph_type[type] == "UW")
        GoTest.template operator()<int>("UW");
    else
    {
        GoTest.template operator()<double>("IW");
        GoTest.template operator()<long double>("IW");
    }       
}

int main(int argc, char* argv[])
{
    Load_Parameters();
    Build_HL(argv[1], argv[2], atoi(argv[3]));
    return 0;
}