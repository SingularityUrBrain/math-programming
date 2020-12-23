#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <fstream>


unsigned int dijkstra(unsigned int st, unsigned int end, const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &graph){
    std::set<std::pair<unsigned int, unsigned int>> cands {std::make_pair(0, st)};
    std::vector<unsigned int> ds(graph.size(), UINT32_MAX);
    ds[st] = 0;
    while (cands.size()){
        auto it = cands.begin();
        unsigned int v = (*it).second;
        for(auto &p: graph[v]){
            if(ds[v] + p.first < ds[p.second]){
                if (ds[p.second] != UINT32_MAX)
                    cands.erase(std::make_pair(ds[p.second], p.second));  // to update
                ds[p.second] = ds[v] + p.first;
                cands.insert(std::make_pair(ds[p.second], p.second));
            }
        }
        cands.erase(it);
    }
    return ds[end];   
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t m, n;
    std::ifstream input("input.txt");
    input >> n >> m;
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> graph(n);
    unsigned int u, v, w;
    for(int i=0;i < m; i++){
        input >> u >> v >> w;
        graph[u-1].push_back(std::make_pair(w, v-1));
        graph[v-1].push_back(std::make_pair(w, u-1));
    }
    input.close();

    unsigned int min_cost = dijkstra(0, n-1, graph);
    std::ofstream output("output.txt");
    output << min_cost;
    output.close();

    return 0;
}