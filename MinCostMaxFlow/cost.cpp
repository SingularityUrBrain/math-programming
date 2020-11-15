#include <iostream>
#include <fstream>
#include <vector>


typedef long long ll;

const int MAXN = 605;
const ll INF = (ll)(1e15);

struct edge {
    int from, to;
    ll flow, cap;
    ll cost;
};



std::vector <ll> net[MAXN];
std::vector <edge> edges;
ll d[MAXN];
bool used[MAXN];
size_t parent[MAXN];
ll phi[MAXN];



void ford_bellman(int n) {
    for (int i = 0; i < n; i++) 
        d[i] = INF;
    d[0] = 0;
    bool ex = false;
    while (!ex) {
        ex = true;
        for (auto &e: edges) {
            int u = e.from, v = e.to;
            if (d[v] > d[u] + e.cost && e.flow < e.cap) {
                d[v] = d[u] + e.cost;
                ex = false;
            }
        }
    }
}

ll astar(int n) {
    for (int i = 0; i < n; i++) {
        d[i] = INF;
        used[i] = false;
    }
    d[0] = 0;
    for (int i = 0; i < n; i++) {
        int u = -1;
        for (int j = 0; j < n; j++)
            if (!used[j] && (u == -1 || d[j] < d[u]))
                u = j;
        used[u] = true;
        for (auto &k: net[u]) {
            int to = edges[k].to;
            ll w = edges[k].cost + phi[u] - phi[to];
            if (d[u] + w < d[to] && edges[k].flow < edges[k].cap) {
                d[to] = d[u] + w;
                parent[to] = k;
            }
        }
    }
    return d[n-1];
}

ll min_cost(int flow, int n) {
    ll min_cost = 0;
    ford_bellman(n);  
    for (int i = 0; i < n; ++i) phi[i] = d[i];
    while (true) {
        if (astar(n) == INF)
            return min_cost;
        for (int i = 0; i < n; ++i)
            phi[i] += d[i];
        ll sent = flow;
        int u = n-1;
        while (u != 0) {
            edge tmp = edges[parent[u]];
            sent = std::min(sent, (tmp.cap - tmp.flow));
            u = tmp.from;
        }
        flow -= sent;
        u = n-1;
        while (u != 0) {
            int v = parent[u];
            edges[v].flow += sent;
            edges[v ^ 1].flow -= sent;
            min_cost += sent*edges[v].cost;
            u = edges[v].from;
        }
        if (flow == 0) break;
    }
    return min_cost;
}   

void add_edge(int from, int to, ll cap, ll cost) {
    net[from].push_back(edges.size());
    edges.push_back({from, to, 0, cap, cost});
    net[to].push_back(edges.size());
    edges.push_back({to, from, cap, cap, -cost});
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cout.tie(nullptr);
    std::fstream cin("input.txt", std::fstream::in);
    std::fstream cout("output.txt", std::fstream::out);

    int n, m, N;
    cin >> n >> m;
    N = n + m + 2;
    std::vector<int> a(n), b(m);
    std::vector<std::vector<ll>> cost (n+1, std::vector<ll>(m+1));

    for(int i = 0; i < n; i++) cin >> a[i];
    for(int i = 0; i < m; i++) cin >> b[i];
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            cin >> cost[i][j];

    // Build net
    for (int i = 1; i <= n; i++) 
        add_edge(0, i, a[i-1], 0);
    for (int i = 1; i <= n; i++)
        for (int j = n+1; j < N-1; j++)
            add_edge(i, j, INF, cost[i][j-n]);
    for (int i = n + 1; i < N-1; i++) 
        add_edge(i, N-1, b[i-n-1], 0);
    int flow = 0;
    for (auto &el: a) flow += el;
    cout << min_cost(flow, N) << '\n';

    return 0;
}