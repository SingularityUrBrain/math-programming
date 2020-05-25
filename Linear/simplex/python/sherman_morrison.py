import operator


def inner_product(v1, v2):
    return sum(map(operator.mul, v1, v2))


def transpose(m):
    return [r for r in zip(*m)]


def matrix_product(m1, m2):
    if len(m1[0]) != len(m2):
        raise ValueError('m1, m2 are not suitable')
    m2_t = transpose(m2)
    return [[inner_product(row1, col2) for col2 in m2_t] for row1 in m1]


def vv_operation(v1, v2, op):
    return list(map(op, v1, v2))


def mm_operation(m1, m2, op):
    '''Element wise action between two matrixes'''
    for i in range(len(m1)):
        for j in range(len(m1)):
            m1[i][j] = op(m1[i][j], m2[i][j])
    return m1


def mc_operation(m, k, op):
    '''Element wise action between matrix and constant'''
    for i in range(len(m)):
        for j in range(len(m)):
            m[i][j] = op(m[i][j], k)
    return m


def main():
    n, k = map(int, input().split())
    k -= 1
    A = []
    B = []
    for _ in range(n):
        A.append(list(map(float, input().split())))
    for _ in range(n):
        B.append(list(map(float, input().split())))
    x_vec = list(map(float, input().split()))

    v_vec = [0]*n
    v_vec[k] = 1
    u_vec = list(map(operator.sub, x_vec, [el[k] for el in A]))
    denum = 1 + inner_product(B[k], u_vec)
    if denum == 0:
        print('NO')
        return
    else:
        print('YES')

    bc_vec = [inner_product(row_b, u_vec) for row_b in B]
    numer = [[bc_el*B[k][j] for j in range(len(B))] for bc_el in bc_vec]
    B_inv = mm_operation(B, mc_operation(
        numer, denum, operator.truediv), operator.sub)

    for row in B_inv:
        print(*row)


if __name__ == "__main__":
    main()
