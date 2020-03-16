import operator
import math
from sys import stdin, stdout


def inner_product(v1, v2):
    return sum(map(operator.mul, v1, v2))


def matrix_vector_product(m, v):
    return [inner_product(row, v) for row in m]


def vector_matrix_product(v, m):
    return [inner_product(v, row) for row in zip(*m)]


def sherman_morrison(B, x, k):
    n = len(B)
    l_vec = matrix_vector_product(B, x)
    l_k = l_vec[k]
    l_vec[k] = -1
    for i in range(n):
        l_vec[i] /= -l_k
    B_inv = [[B[j][i] + l_vec[j]*B[k][i] if j != k else l_vec[j]*B[k][i]
              for i in range(n)] for j in range(n)]
    return B_inv


def invert_matrix(A):
    n = len(A)
    E = [[0]*n for i in range(n)]
    #A = copy.deepcopy(A)
    for i in range(n):
        E[i][i] = 1
    for i in range(n):
        maxrow = i
        for j in range(i+1, n):
            if abs(A[j][i]) > abs(A[maxrow][i]):
                maxrow = j
        if maxrow != i:
            A[i], A[maxrow] = A[maxrow], A[i]
            E[i], E[maxrow] = E[maxrow], E[i]
        for j in range(i+1, n):
            c = A[j][i]/A[i][i]
            for k in range(i):
                E[j][k] -= E[i][k]*c
            for k in range(i, n):
                A[j][k] -= A[i][k]*c
                E[j][k] -= E[i][k]*c
    for i in range(n-1, -1, -1):
        c = A[i][i]
        for j in range(i):
            for k in range(n-1, -1, -1):
                E[j][k] -= E[i][k]*A[j][i]/c
        for j in range(n):
            E[i][j] /= c
    return E


def main():
    inp = stdin.readline
    m, n = map(int, inp().split())
    A = [tuple(map(float, inp().split())) for _ in range(m)]
    inp()  # b should be here
    c = tuple(map(float, inp().split()))
    x = list(map(float, inp().split()))
    j_basis = list(map(int, inp().split()))

    c_basis = [c[i-1] for i in j_basis]
    A_basis = [[A[i][j-1] for j in j_basis] for i in range(m)]
    j_n = [i for i in range(n) if i+1 not in j_basis]
    B = invert_matrix(A_basis)
    while True:
        u_vec = vector_matrix_product(c_basis, B)
        uA_vec = vector_matrix_product(u_vec, A)
        min_delta = math.inf
        for j in j_n:
            d = uA_vec[j] - c[j]
            if d < min_delta:
                min_delta = d
                j0 = j
        if min_delta >= 0:
            stdout.write('Bounded\n')
            stdout.write(' '.join(str(el) for el in x))
            stdout.write('\n')
            return
        col_to_replace = [row[j0] for row in A]
        z_vec = matrix_vector_product(B, col_to_replace)
        omega_min = math.inf
        for j in range(m):
            if z_vec[j] > 0:
                om_i = x[j_basis[j]-1]/z_vec[j]
                if om_i < omega_min:
                    omega_min = om_i
                    s = j
        if omega_min == math.inf:
            stdout.write('Unbounded\n')
            return
        # new plan
        for j in j_n:
            x[j] = 0
        x[j0] = omega_min
        for ind, j in enumerate(j_basis):
            x[j-1] -= omega_min*z_vec[ind]

        j_n.remove(j0)
        j_n.append(j_basis[s]-1)
        # new basis
        j_basis[s] = j0 + 1
        c_basis[s] = c[j0]

        B = sherman_morrison(B, col_to_replace, s)


if __name__ == "__main__":
    main()
