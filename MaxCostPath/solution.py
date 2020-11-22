from collections import deque
import sys

def maxcost(graph, src, t):
    used = set()
    ldag = deque()
    def topological_sort(u):
        used.add(u)
        for v, c in graph[u]:
            if v not in used:
                topological_sort(v)
        ldag.append(u)
    topological_sort(src)
    cost = [-1]*len(graph)
    cost[src] = 0
    while ldag:
        u = ldag.pop()
        if u == t:
            return cost[u]
        if cost[u] != -1:
            for v, c in graph[u]:
                d = cost[u] + c
                if cost[v] < d:
                    cost[v] = d
    return cost[t]

 
if __name__ == '__main__':
    inp = sys.stdin.readline
    n, m = map(int, inp().split())
    graph = [[] for _ in range(n)]
    for i in range(m):
        u, v, c = map(int, inp().split())
        graph[u-1].append((v-1, c))
    sys.stdout.write(str(maxcost(graph, 0, n-1)))
 