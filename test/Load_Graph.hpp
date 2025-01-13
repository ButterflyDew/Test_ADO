#ifndef LOADGRAPH_H
#define LOADGRAPH_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "../AdjacencyList.hpp"
template<typename NodeType, typename DistanceType>
class LoadGraph
{
    public:
        int N;
        AdjacencyList<NodeType, DistanceType> adjList;
        LoadGraph(std::string filepath, int isweight, int offset)
        {
            std::string name[2] = {"/graph.txt","/Weightgraph.txt"};
            std::string filename = filepath + name[isweight];
            
            std::ifstream inputFile(filename);

            if (!inputFile.is_open()) 
                std::cerr << "Can't open file: " << filename << std::endl;

            std::string line;
            getline(inputFile, line);

            N = extractIntegers(line)[0];
            adjList = AdjacencyList<NodeType, DistanceType>(N);

            while (getline(inputFile, line)) 
            {
                if(isweight == 0)
                {
                    auto info = extractIntegers(line);
                    adjList.AddUndirectedEdge(info[0] - offset, info[1] - offset, 1);
                }
                else
                {
                    std::istringstream iss(line);
                    int u,v; double w;
                    iss>>u>>v>>w;
                    u -= offset, v -= offset;
                    adjList.AddUndirectedEdge(u, v, w);
                }
            }
            inputFile.close();
        }
    
    private:
        std::vector<int> extractIntegers(const std::string& input) 
        {
            std::vector<int> integers;
            std::istringstream iss(input);
            std::string token;

            while (iss >> token) 
            {
                try 
                {
                    int number = stoi(token);
                    integers.push_back(number);
                } catch (const std::invalid_argument& e) 
                {
                    std::cerr << "not a number" << std::endl;
                }
            }

            return integers;
        }
};

#endif