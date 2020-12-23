from collections import deque

INF = float('inf')


def bfs(net, flow, s, t):
    q = deque([s])
    paths = {s: []}
    if s == t:
        return paths[s]
    while q:
        u = q.popleft()
        for v in range(len(net)):
            if v not in paths and (net[u][v] - flow[u][v] > 0):
                paths[v] = paths[u] + [(u, v)]
                if v == t:
                    return paths[v]
                q.append(v)


def main():
    seq = []
    with open('flow.in') as f:
        m, _ = map(int, f.readline().split())
        net = [[0]*m for _ in range(m)]
        flow = [[0]*m for _ in range(m)]
        s, t = map(int, f.readline().split())
        s -= 1
        t -= 1
        for line in f:
            u, v, c = map(int, line.split())
            u -= 1
            v -= 1
            seq.append((u, v, c))
            net[u][v] = c
            
    max_flow = 0
    path = bfs(net, flow, s, t)
    while path:
        sent = min(net[u][v] - flow[u][v] for u, v in path)
        max_flow += sent
        for u, v in path:
            flow[u][v] += sent
            flow[v][u] -= sent
        path = bfs(net, flow, s, t)
    
    with open('flow.out', 'w') as out:
        out.write(str(max_flow))
        for u, v, c in seq:
            out.write(f'\n{max(0, flow[u][v])}')


if __name__ == "__main__":
    main()
