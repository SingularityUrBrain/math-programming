import math


def find_ifs(a, b):
    '''The North-west method to find initial feasible solution'''
    m, n = len(a), len(b)
    X = [[0]*n for _ in range(m)]
    basis = set()
    i = j = 0
    while i < m and j < n:
        basis.add((i, j))
        if a[i] < b[j]:
            X[i][j] = a[i]
            b[j] -= a[i]
            i += 1
        else:
            X[i][j] = b[j]
            a[i] -= b[j]
            j += 1

    no_basis = set((i, j) for j in range(n)
                   for i in range(m) if (i, j) not in basis)
    if len(basis) != m + n - 1:
        for i, j in no_basis:
            basis.add((i, j))
            if find_cycle(m, n, basis) is None:
                if len(basis) == m + n - 1:
                    break
            else:
                basis.remove((i, j))
        no_basis.difference_update(basis)
    return X, basis, no_basis


def find_cycle(m, n, basis):
    basis = basis.copy()
    rows, cols = [0]*m, [0]*n
    for i, j in basis:
        rows[i] += 1
        cols[j] += 1
    while True:
        is_ex = True
        for k in range(m):
            if rows[k] == 1:
                is_ex = False
                for i, j in basis:
                    if i == k:
                        cols[j] -= 1
                        rows[i] = 0
                        basis.remove((i, j))
                        break
        for k in range(n):
            if cols[k] == 1:
                is_ex = False
                for i, j in basis:
                    if j == k:
                        rows[i] -= 1
                        cols[j] = 0
                        basis.remove((i, j))
                        break
        if is_ex:
            return basis
        if len(basis) < 4:
            return None


def find_potentials(C, basis):
    inf = float('inf')
    u = [inf]*len(C)
    v = [inf]*len(C[0])
    u[0] = 0
    for _ in range(len(basis)):
        for i, j in basis:
            if v[j] == inf and u[i] != inf:
                v[j] = C[i][j] - u[i]
                break
            elif u[i] == inf and v[j] != inf:
                u[i] = C[i][j] - v[j]
                break
    return u, v


def split_cycle(i0, j0, cycle):
    neg, pos = set(), set()
    pos.add((i0, j0))
    for _ in range(len(cycle) >> 1):
        for i, j in cycle:
            if i == i0 and j != j0:
                neg_i, neg_j = i, j
                break
        neg.add((neg_i, neg_j))
        for i, j in cycle:
            if j == neg_j and i != neg_i:
                i0, j0 = i, j
                break
        pos.add((i0, j0))
    return neg, pos


def solve(a, b, C):
    m, n = len(a), len(b)
    diff = sum(a) - sum(b)
    if diff < 0:
        a.append(abs(diff))
        m += 1
        C.append([0 for _ in range(n)])
    elif diff > 0:
        b.append(diff)
        n += 1
        for row in C:
            row.append(0)
    X, basis, no_basis = find_ifs(a, b)
    if len(no_basis) == 0:
        return X

    while True:
        u, v = find_potentials(C, basis)
        i0, j0 = min(no_basis, key=lambda x: C[x[0]][x[1]]-u[x[0]]-v[x[1]])
        min_delta = C[i0][j0]-u[i0]-v[j0]
        if min_delta >= 0:
            return X

        basis.add((i0, j0))
        no_basis.remove((i0, j0))
        cycle = find_cycle(m, n, basis)
        neg, pos = split_cycle(i0, j0, cycle)
        i_star, j_star = min(neg, key=lambda el: X[el[0]][el[1]])
        theta = X[i_star][j_star]
        for el in pos:
            X[el[0]][el[1]] += theta
        for el in neg:
            X[el[0]][el[1]] -= theta
        basis.remove((i_star, j_star))
        no_basis.add((i_star, j_star))


def main():
    m, n = map(int, input().split())
    # cost
    C = [list(map(int, input().split())) for _ in range(m)]
    # capacity
    a = list(map(int, input().split()))
    # demand
    b = list(map(int, input().split()))
    X = solve(a, b, C)
    for row in X[:m]:
        print(*row[:n])


if __name__ == "__main__":
    main()
