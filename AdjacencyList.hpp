#ifndef ORACLE_H
#define ORACLE_H

#include <limits>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <unordered_set>

template<typename NodeType, typename DistanceType>
class AdjacencyList
{
    public:
        using Edges = std::vector<std::pair<NodeType, DistanceType>>;
        using Matrix = std::vector<Edges>;

        AdjacencyList(){}
        /*Init |V|, 0-index*/
        AdjacencyList(std::size_t v) : sz(v), edgeCount(0) { mat.assign(v, std::vector<std::pair<NodeType, DistanceType>>()); }

        void AddUndirectedEdge(NodeType u, NodeType v, DistanceType d)
        {
            mat[u].push_back(std::make_pair(v, d));
            mat[v].push_back(std::make_pair(u, d));

            auto uv = std::make_pair(u, v), vu = std::make_pair(v, u);
            if(weig.contains(uv)) weig[uv] = std::min(weig[uv], d);
            else weig[uv] = d;
            if(weig.contains(vu)) weig[vu] = std::min(weig[vu], d);
            else weig[vu] = d;

            edgeCount++;
        }

        const std::vector<std::pair<NodeType, DistanceType>>& GetAllEdges(NodeType v) { return mat[v]; }
        const Matrix& GetMatrix() const { return mat; }
        std::size_t GetSize() const { return sz; }
        std::size_t GetEdgeCount() const { return edgeCount; }
        DistanceType GetEdgeuv(NodeType u, NodeType v)
        {
            auto uv = std::make_pair(u, v);
            if(weig.contains(uv)) return weig[uv];
            else return std::numeric_limits<DistanceType>::max(); 
        }

        // Use some shortest path algorithm
        std::vector<DistanceType> GetNearest(NodeType u) const
        {
            //std::unordered_set<NodeType> vs(v.begin(), v.end());

            std::vector<DistanceType> dist(sz, std::numeric_limits<DistanceType>::max());
            dist[u] = 0;

            std::priority_queue<std::pair<DistanceType, NodeType>, std::vector<std::pair<DistanceType, NodeType>>, std::greater<std::pair<DistanceType, NodeType>>> q;
            q.push({0, u});

            while(!q.empty())
            {
				auto temp = q.top();
				auto d = temp.first;auto a =temp.second;
                q.pop();

                if(dist[a] < d) continue;

 				for(auto temp1: mat[a])
                {
					auto b = temp1.first;auto d1 =temp1.second;
                    if(d1 + d < dist[b])
                    {
                        dist[b] = d1 + d;
                        q.push({d1 + d, b});
                    }
                }
            }
            return dist;
        }

		std::vector<DistanceType> GetDistance(NodeType u) const
        {
            return GetNearest(u);
        }

    private:
        std::size_t sz;
        std::size_t edgeCount;
        Matrix mat;
        std::map <std::pair <NodeType, NodeType>, DistanceType> weig; 
};


#endif