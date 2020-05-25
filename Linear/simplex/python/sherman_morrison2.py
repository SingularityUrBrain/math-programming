import operator


def inner_product(v1, v2):
    return sum(map(operator.mul, v1, v2))


def matrix_vector_product(m, v):
    return [inner_product(row, v) for row in m]


def sherman_morrison(B, x, k):
    n = len(B)
    l_vec = matrix_vector_product(B, x)
    l_k = l_vec[k]
    if l_k == 0:
        # raise ValueError('Matrix is singular.')
        return
    l_vec[k] = -1
    l_vec = [-el/l_k for el in l_vec]
    B_inv = [[B[j][i] + l_vec[j]*B[k][i] if j != k else l_vec[j]*B[k][i]
              for i in range(n)] for j in range(n)]
    return B_inv


def main():
    n, k = map(int, input().split())
    for _ in range(n):
        input()
    B = [list(map(float, input().split())) for _ in range(n)]
    x_vec = list(map(float, input().split()))
    B_inv = sherman_morrison(B, x_vec, k-1)
    if B_inv:
        print('YES')
        for row in B_inv:
            print(*row)
    else:
        print('NO')


if __name__ == "__main__":
    main()
