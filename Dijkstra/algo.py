import heapq as pq

## O(n + m*log(m))
def dijkstra(st, end, graph):
    heap = [(0, st)]
    ds = [float('inf')]*len(graph)
    ds[st] = 0
    used = set()
    while heap:
        d, v = pq.heappop(heap)
        if d <= ds[v] and v not in used: 
            used.add(v)
            for to_len, to in graph[v]:
                if ds[v] + to_len < ds[to]:
                    ds[to] = ds[v] + to_len
                    pq.heappush(heap, (ds[to], to))
    return ds[end]


if __name__ == "__main__":    
    with open('input.txt') as f:
        n, m = map(int, f.readline().split())
        graph = [[] for i in range(n)]
        for line in f:
            u, v, w = map(int, line.split())
            graph[u-1].append((w, v-1))
            graph[v-1].append((w, u-1))
    with open('output.txt', 'w') as f:
        f.write(str(dijkstra(0, n-1, graph)))
