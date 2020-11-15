from collections import deque

INF = float('inf')

class Edge:
    def __init__(self, v, reverse, flow, cap):
        self.v = v
        self.reverse = reverse
        self.flow = flow
        self.cap = cap

class Graph:
    def __init__(self, n):
        self.adj = [[] for _ in range(n)]
        self.d = [0]*n
        self.n = n
        self.edges = []

    def add_edge(self, u, v, c):
        self.adj[u].append(len(self.edges))
        self.adj[v].append(len(self.edges)+1)
        self.edges.append(Edge(v, u, 0, c))
        self.edges.append(Edge(u, v, 0, 0))

    
    def dfs(self, u, t, fmin):
        if u == t:
            return fmin              
        for i in self.adj[u]:
            edge = self.edges[i]
            capacity = edge.cap - edge.flow
            if self.d[edge.v] == self.d[u] + 1 and capacity > 0:
                delta = self.dfs(edge.v, t, min(fmin, capacity))
                if delta > 0:
                    edge.flow += delta
                    self.edges[i^1].flow -= delta
                    return delta
        return 0

    def bfs(self, s, t):
        for i in range(self.n): 
            self.d[i] = INF
        q = deque([s])
        self.d[s] = 0
        while q:
            u = q.popleft()
            for i in self.adj[u]:
                edge = self.edges[i]
                if self.d[edge.v] == INF and edge.flow < edge.cap:
                    self.d[edge.v] = self.d[u] + 1
                    q.append(edge.v)
        return self.d[t] != INF

    def get_max_flow(self, s, t):
        max_flow = 0
        while self.bfs(s, t):
            while True:
                sent = self.dfs(s, t, INF)
                if sent == 0:
                    break
                max_flow += sent
        return max_flow


def main():
    with open('flow.in') as f:
        m, n = map(int, f.readline().split())
        s, t = map(int, f.readline().split())
        net = Graph(m)
        for line in f:
            u, v, c = map(int, line.split())
            net.add_edge(u-1, v-1, c)
    
    max_flow = net.get_max_flow(s-1, t-1)
    with open('flow.out', 'w') as out:
        out.write(str(max_flow))
        for e in net.edges[::2]:
            out.write(f'\n{e.flow}')
    

if __name__ == "__main__":
    main()
