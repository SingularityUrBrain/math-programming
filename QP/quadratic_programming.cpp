#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <cmath>

#define INF 1e10
#define EPS 1e-5

template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>> &m, const std::vector<T> &v)
{
    std::vector<T> product_vec(m.size());
    for (size_t i = 0; i < m.size(); i++)
    {
        for (size_t j = 0; j < v.size(); j++)
            product_vec[i] += m[i][j] * v[j];
    }
    return product_vec;
}
template <typename T>
std::vector<T> operator*(const std::vector<T> &v, const std::vector<std::vector<T>> &m)
{
    std::vector<T> product_vec(m[0].size());
    for (size_t i = 0; i < m[0].size(); i++)
    {
        for (size_t j = 0; j < v.size(); j++)
            product_vec[i] += v[j] * m[j][i];
    }
    return product_vec;
}

template <typename T>
std::vector<T> operator+(const std::vector<T> &v1, const std::vector<T> &v2)
{
    std::vector<T> res(v1.size());
    for (int i = 0; i < v1.size(); i++)
        res[i] = v1[i] + v2[i];
    return res;
}
template <typename T>
std::vector<T> operator-(const std::vector<T> &v1, const std::vector<T> &v2)
{
    std::vector<T> res(v1.size());
    for (int i = 0; i < v1.size(); i++)
        res[i] = v1[i] - v2[i];
    return res;
}
template <typename T>
T operator*(const std::vector<T> &v1, const std::vector<T> &v2)
{
    T res = T();
    for (size_t i = 0; i < v1.size(); i++)
        res += v1[i] * v2[i];
    return res;
}
template <typename T>
std::vector<T> operator-(const std::vector<T> &v)
{
    std::vector<T> neg_vec(v.size());
    for (size_t i = 0; i < v.size(); i++)
        neg_vec[i] = -v[i];
    return neg_vec;
}
template <typename T>
std::vector<T> operator*(T num, const std::vector<T> &v)
{
    std::vector<T> res(v.size());
    for (size_t i = 0; i < v.size(); i++)
        res[i] = num * v[i];
    return res;
}
template <typename T>
std::vector<T> operator*(const std::vector<T> &v, T num)
{
    std::vector<T> res(v.size());
    for (size_t i = 0; i < v.size(); i++)
        res[i] = num * v[i];
    return res;
}

template <typename T>
std::vector<std::vector<T>> invert_matrix(std::vector<std::vector<T>> &A)
{
    std::vector<std::vector<T>> E(A.size(), std::vector<T>(A.size()));
    for (size_t i = 0; i < E.size(); i++)
        ++E[i][i];

    for (size_t i = 0; i < A.size(); i++)
    {
        size_t max_row = i;
        for (size_t j = i + 1; j < A.size(); j++)
            if (std::fabs(A[j][i]) > std::fabs(A[max_row][i]))
                max_row = j;
        if (max_row != i)
        {
            A[i].swap(A[max_row]);
            E[i].swap(E[max_row]);
        }
        for (size_t j = i + 1; j < A.size(); j++)
        {
            T c = A[j][i] / A[i][i];
            for (size_t k = 0; k < i; k++)
                E[j][k] -= E[i][k] * c;
            for (size_t k = i; k < A.size(); k++)
            {
                A[j][k] -= A[i][k] * c;
                E[j][k] -= E[i][k] * c;
            }
        }
    }
    for (int i = A.size() - 1; i >= 0; i--)
    {
        T c = A[i][i];
        for (size_t j = 0; j < i; j++)
            for (int k = E.size() - 1; k >= 0; k--)
                E[j][k] -= E[i][k] * A[j][i] / c;
        for (size_t j = 0; j < E.size(); j++)
            E[i][j] /= c;
    }
    return E;
}

template <typename T>
std::vector<T> solve_qp(const std::vector<std::vector<T>> &A, const std::vector<T> &c,
                        const std::vector<std::vector<T>> &D, std::vector<T> x, std::set<size_t> j_sb, std::set<size_t> j_eb)
{
    size_t m = A.size(), n = A[0].size(), k = j_eb.size();
    std::vector<std::vector<T>> A_sb(m, std::vector<T>(m)), B;
    std::vector<T> c_r, c_sb(m), u, l_dir(n);
    std::set<size_t> j_neb;

    for (size_t i = 0; i < n; i++)
    {
        if (j_eb.count(i) == 0)
            j_neb.insert(i);
    }
    auto it_set = j_sb.begin();
    for (size_t i = 0; i < m; i++, it_set++)
        for (size_t j = 0; j < m; j++)
            A_sb[j][i] = A[j][*it_set];
    B = invert_matrix<T>(A_sb);
    bool is_need_change;

    while (true)
    {
        is_need_change = false;
        c_r = D * x + c;
        it_set = j_sb.begin();
        for (size_t i = 0; i < m; i++, it_set++)
        {
            c_sb[i] = c_r[*it_set];
        }
        u = -c_sb * B;
        T min_delta = INF;
        size_t j0;
        for (auto &&it : j_neb)
        {
            T pr = T();
            for (size_t j = 0; j < m; j++)
            {
                pr += u[j] * A[j][it];
            }
            pr += c_r[it];
            if (pr < -EPS)
            {
                min_delta = pr;
                j0 = it;
                break;
            }
        }
        if (min_delta >= 0)
            return x;
        // step -- 3
        while (true)
        {
            for (auto &&j : j_neb)
                l_dir[j] = 0;
            l_dir[j0] = 1.;
            // fill H
            std::vector<std::vector<T>> H(m + k, std::vector<T>(m + k));
            it_set = j_eb.begin();
            for (size_t i = 0; i < k; i++, it_set++)
            {
                auto it_set2 = j_eb.begin();
                for (size_t j = 0; j < k; j++, it_set2++)
                {
                    H[i][j] = D[*it_set][*it_set2];
                }
            }
            it_set = j_eb.begin();
            for (size_t i = 0; i < k; i++, it_set++)
            {
                for (size_t j = k, a_i = 0; j < k + m; j++, a_i++)
                {
                    H[i][j] = A[a_i][*it_set];
                }
            }
            for (size_t i = k, a_i = 0; i < k + m; i++, a_i++)
            {
                it_set = j_eb.begin();
                for (size_t j = 0; j < k; j++, it_set++)
                {
                    H[i][j] = A[a_i][*it_set];
                }
            }

            // fill bb
            std::vector<T> bb;
            for (auto &&it : j_eb)
            {
                bb.push_back(D[it][j0]);
            }
            for (size_t i = 0; i < m; i++)
            {
                bb.push_back(A[i][j0]);
            }
            // fill l*
            std::vector<T> Hbb = -invert_matrix(H) * bb;
            it_set = j_eb.begin();
            for (size_t i = 0; i < k; i++, it_set++)
                l_dir[*it_set] = Hbb[i];

            // Find min Omega
            T min_theta = INF;
            size_t j_star;
            for (auto &&it : j_eb)
            {
                if (l_dir[it] < -EPS)
                {
                    T el = -x[it] / l_dir[it];
                    if (el < min_theta)
                    {
                        min_theta = el;
                        j_star = it;
                    }
                }
            }
            T d = bb * Hbb + D[j0][j0];
            if (d > EPS)
            {
                T el = std::fabs(min_delta) / d;
                if (el < min_theta)
                {
                    min_theta = el;
                    j_star = j0;
                }
            }
            if (min_theta == INF)
            {
                std::cout << "Unbounded\n";
                exit(0);
            }
            // update x
            x = x + min_theta * l_dir;

            // case 1
            if (j_star == j0)
            {
                j_eb.insert(j0);
                j_neb.erase(j0);
                ++k;
                if(is_need_change)
                {
                    auto it_set = j_sb.begin();
                    for (size_t i = 0; i < m; i++, it_set++)
                        for (size_t j = 0; j < m; j++)
                            A_sb[j][i] = A[j][*it_set];
                    B = invert_matrix<T>(A_sb);
                }
                break;
            }
            else
            {
                std::set<size_t> j_e_nosb;
                std::set_difference(j_eb.begin(), j_eb.end(), j_sb.begin(), j_sb.end(),
                                    std::inserter(j_e_nosb, j_e_nosb.end()));
                // case 2
                if (j_e_nosb.count(j_star) == 1)
                {
                    j_eb.erase(j_star);
                    j_neb.insert(j_star);
                    --k;
                    min_delta += min_theta * d;
                }
                // case 3,4
                else
                {
                    bool is_case3 = false;
                    size_t j_plus;
                    size_t s = std::distance(j_sb.begin(), j_sb.find(j_star));
                    for (auto &&j_p : j_e_nosb)
                    {
                        T pr = T();
                        for (size_t i = 0; i < m; i++)
                        {
                            pr += B[s][i] * A[i][j_p];
                        }
                        if (std::fabs(pr) > EPS)
                        {
                            is_case3 = true;
                            j_plus = j_p;
                            break;
                        }
                    }

                    j_sb.erase(j_star);
                    j_eb.erase(j_star);
                    j_neb.insert(j_star);
                    // case 3
                    if (is_case3)
                    {
                        j_sb.insert(j_plus);
                        --k;
                        min_delta += min_theta * d;
                        is_need_change = true;
                    }
                    // case 4
                    else
                    {
                        j_sb.insert(j0);
                        j_eb.insert(j0);
                        j_neb.erase(j0);
                        auto it_set = j_sb.begin();
                        for (size_t i = 0; i < m; i++, it_set++)
                            for (size_t j = 0; j < m; j++)
                                A_sb[j][i] = A[j][*it_set];
                        B = invert_matrix<T>(A_sb);
                        break;
                    }
                }
            }
        }
    }
    return x;
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t m, n, k;
    std::cin >> m >> n >> k;
    std::vector<std::vector<double>> A(m, std::vector<double>(n)), D(n, std::vector<double>(n));
    std::vector<double> c(n), b(m), x(n);
    std::set<size_t> j_sb, j_eb;
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
            std::cin >> A[i][j];
    }
    for (size_t i = 0; i < m; i++)
        std::cin >> b[i];
    for (size_t i = 0; i < n; i++)
        std::cin >> c[i];
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
            std::cin >> D[i][j];
    }
    for (size_t i = 0; i < n; i++)
        std::cin >> x[i];
    size_t index;
    for (size_t i = 0; i < m; i++)
    {
        std::cin >> index;
        j_sb.insert(index - 1);
    }
    for (size_t i = 0; i < k; i++)
    {
        std::cin >> index;
        j_eb.insert(index - 1);
    }

    auto x0 = solve_qp<double>(A, c, D, x, j_sb, j_eb);
    std::cout << "Bounded" << '\n';
    std::cout << std::setprecision(10) << std::fixed;
    for (auto &&el : x0)
        std::cout << el << ' ';
    std::cout << '\n';
    return 0;
}

// 3 8 3
// 1 2 0 1 0 4 -1 -3
// 1 3 0 0 1 -1 -1 2
// 1 4 1 0 0 2 -2 0
// 4 5 6
// -10 -31 7 0 -21 -16 11 -7
// 6 11 -1 0 6 -7 -3 -2
// 11 41 -1 0 7 -24 0 -3
// -1 -1 1 0 -3 -4 2 -1
// 0 0 0 0 0 0 0 0
// 6 7 -3 0 11 6 -7 1
// -7 -24 -4 0 6 42 -7 10
// -3 0 2 0 -7 -7 5 -1
// -2 -3 -1 0 1 10 -1 3
// 0 0 6 4 5 0 0 0
// 3 4 5
// 3 4 5
