#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#define INF 1e9

int min_dist(std::vector<std::vector<std::vector<int>>> &dp, std::vector<std::vector<int>> &dist, int k) {

    int mask = 1 << dist.size();
    for(int s = 1; s < mask; s++)
    {
        for(int i = 0; i < dist.size(); i++){
            if (s & (1<<i)){
                for(int j = 0; j < dist.size(); j++){
                    if (((1<<j) & s) && i != j){
                        int prev_state = s ^ (1 << j);
                        for(int t = 0; t <= k; t++){
                            dp[s][j][t] = std::min(dp[s][j][t], dist[i][j] + dp[prev_state][i][t]);
                            if(t < k)
                                dp[s][j][t+1] = std::min(dp[s][j][t+1], dp[prev_state][i][t]);
                        }
                    } 
                }
            }
        }
    }
    int ans = INT32_MAX;
    for (int i = 0; i < dist.size(); i++) 
        ans = std::min(ans, dp[mask-1][i][k]);
    return ans >= INF ? -1 : ans;
}


int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t m, n, k;
    std::cin >> n >> m >> k;
    std::vector<std::vector<int>> graph(n, std::vector<int>(n, INF));
    std::map<std::string, int> map;
    for (int i = 0; i < n; i++){
        std::string x;
        std::cin >> x;
        map.insert(std::pair<std::string, int>(x, i));
    }
    std::string c1, c2;
    int d;
    for (size_t i = 0; i < m; i++){
        std::cin >> c1 >> c2 >> d; 
        graph[map[c1]][map[c2]] = d;
        graph[map[c2]][map[c1]] = d;
    }
    std::vector<std::vector<std::vector<int>>> dp(1<<n, std::vector<std::vector<int>>(n, std::vector<int>(k+1, INF)));
    // base case
    for (int i = 0; i < n; i++){
        dp[1<<i][i][0] = 0;
        for (int j = 0; j <= k; j++)
            dp[1<<i][i][j] = 0;
    }

    std::cout << min_dist(dp, graph, k) << '\n';
    
    return 0;
}