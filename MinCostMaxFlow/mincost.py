INF = float('inf')


def ford_bellman(net, cost, d, u, N):
    for i in range(N):
        d[i] = INF
    d[u] = 0
    ex = False
    while not ex:
        ex = True
        for u in range(N):
            for v, cap in enumerate(net[u]):
                if d[v] > d[u] + cost[u][v] and cap > 0:
                    d[v] = d[u] + cost[u][v]
                    ex = False

def astar(net, cost, d, phi, used, parent, s, N):
    for i in range(N):
        d[i] = INF
        used[i] = False
    d[s] = 0
    for i in range(N):
        u = -1
        for j in range(N): 
            if not used[j] and (d[j] - d[u] < 0 or u == -1):
                u = j
        used[u] = True
        for v, cap in enumerate(net[u]):
            w = cost[u][v] + phi[u] - phi[v]
            if d[u] + w < d[v] and cap > 0:
                d[v] = d[u] + w
                parent[v] = u



def main():
    with open('input.txt') as f:
        n, m = map(int, f.readline().split())
        N = n + m + 2
        s = 0
        t = N - 1
        a = list(map(int, f.readline().split()))
        flow = sum(a)
        b = list(map(int, f.readline().split()))
        net = [[0]*N for _ in range(N)]
        
        # Build net
        for i in range(1, n+1):
            net[0][i] = a[i-1]
        for i in range(1, n+1):
            for j in range(1 + n, N - 1):
                net[i][j] = float('inf')
        j = 0
        for i in range(1 + n, N - 1):
            net[i][N-1] = b[j]
            j += 1

        # Build cost matrix
        cost = [[0]*N]
        for line in f:
            c = [0]*(n+1) + list(map(int, line.split())) + [0]
            cost.append(c)
        cost.extend([[0]*N for _ in range(m+1)])
        for i in range(1 + n, N - 1):
            for j in range(1, n + 1):
                cost[i][j] = -cost[j][i]

    parent = [INF]*N
    d = [INF]*N
    seen = [False]*N
    min_cost = 0

    ford_bellman(net, cost, d, s, N)
    phi = d[::]
    while True:
        astar(net, cost, d, phi, seen, parent, s, N)
        if d[t] == INF:
            break
        for i in range(N):
            phi[i] += d[i]
        sent = flow
        v = t
        while v != s:
            sent = min(sent, net[parent[v]][v])
            v = parent[v]
        flow -= sent
        v = t
        while v != s:
            u = parent[v]
            net[u][v] -= sent
            net[v][u] += sent
            min_cost += sent*cost[u][v]
            v = u
        if flow == 0:
            break
    with open('output.txt', 'w') as out:
        out.write(str(min_cost))


if __name__ == "__main__":
    main()

