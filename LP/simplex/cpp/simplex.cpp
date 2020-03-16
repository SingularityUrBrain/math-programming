#include <iostream> 
#include <iomanip>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <queue>

#define INF 1e16
#define EPS 1e-8

template<typename T>
std::vector<T> matrix_vector_product(const std::vector<std::vector<T>>& m, const std::vector<T>& v){
    std::vector<T> product_vec(m.size());
    for (size_t i = 0; i < m.size(); i++){
        for (size_t j = 0; j < v.size(); j++)
            product_vec[i] += m[i][j]*v[j];
    }
    return product_vec;
}

template<typename T>
std::vector<T> vector_matrix_product(const std::vector<T>& v, const std::vector<std::vector<T>>& m){
    std::vector<T> product_vec(m[0].size());
    for (size_t i = 0; i < m[0].size(); i++){
        for (size_t j = 0; j < v.size(); j++)
            product_vec[i] += v[j]*m[j][i];
    }
    return product_vec;
}

template<typename T>
std::vector<std::vector<T>> sherman_morrison(const std::vector<std::vector<T>>& B, const std::vector<T>& x, size_t k){
    std::vector<T> l_vec = matrix_vector_product<T>(B, x);
    T l_k = l_vec[k];
    l_vec[k] = -1;
    for (auto &el : l_vec) el /= -l_k;
    std::vector<std::vector<T>> B_inv(B.size(), std::vector<T>(B.size()));
    for (size_t i = 0; i < B.size(); i++)
        if (i != k){
            for (size_t j = 0; j < B.size(); j++)
                B_inv[i][j] = B[i][j] + l_vec[i]*B[k][j];
        }
        else{
            for (size_t j = 0; j < B.size(); j++)
                B_inv[i][j] = l_vec[k]*B[k][j];
        }
    return B_inv;
}

template<typename T>
std::vector<std::vector<T>> invert_matrix(std::vector<std::vector<T>>& A){
    std::vector<std::vector<T>> E(A.size(), std::vector<T>(A.size()));
    for (size_t i = 0; i < E.size(); i++) ++E[i][i];

    for (size_t i = 0; i < A.size(); i++){
        size_t max_row = i;
        for (size_t j = i+1; j < A.size(); j++)
            if (abs(A[j][i]) > abs(A[max_row][i]))
                max_row = j;
        if (max_row != i){
            A[i].swap(A[max_row]);
            E[i].swap(E[max_row]);
        }
        for(size_t j = i+1; j < A.size(); j++){
            T c = A[j][i]/A[i][i];
            for(size_t k = 0; k < i; k++)
                E[j][k] -= E[i][k]*c;
            for(size_t k = i; k < A.size(); k++){
                A[j][k] -= A[i][k]*c;
                E[j][k] -= E[i][k]*c;
            }
        }
    }
    for(int i = A.size()-1; i >= 0; i--){
        T c = A[i][i];
        for(size_t j = 0; j < i; j++)
            for(int k = E.size()-1; k >= 0; k--)
                E[j][k] -= E[i][k]*A[j][i]/c;
        for(size_t j = 0; j < E.size(); j++)
            E[i][j] /= c;
    }
    return E;
}

template<typename T>
std::pair<std::vector<T>, std::vector<size_t>> main_phase(const std::vector<std::vector<T>>& A,
 const std::vector<T>& c, std::vector<T> x, std::vector<size_t> j_basis)
{
    size_t m = A.size(), n = A[0].size();
    std::vector<std::vector<T>> A_basis(m, std::vector<T>(m)), B;
    std::vector<T> c_basis(m);
    std::vector<size_t> j_n;

    for (size_t i = 0; i < m; i++)
        c_basis[i] = c[j_basis[i]-1];
    for (size_t i = 0; i < m; i++)
        for(size_t j = 0; j < m; j++)
            A_basis[j][i] = A[j][j_basis[i]-1];

    std::unordered_set<size_t> j_basis_set(j_basis.begin(), j_basis.end());
    for (size_t i = 0, j = 0; i < n; i++){
        if(!j_basis_set.count(i+1)) j_n.push_back(i);
    }
    B = invert_matrix<T>(A_basis);
    std::vector<T> u_vec, uA_vec, col_to_replace(m), z_vec;

    while (1)
    {
        u_vec = vector_matrix_product<T>(c_basis, B);
        uA_vec = vector_matrix_product<T>(u_vec, A);
        T min_delta = INF;
        size_t j0, j_n_index;
        for (size_t j = 0; j < j_n.size(); j++){
            T d = uA_vec[j_n[j]] - c[j_n[j]];
            if (d < min_delta){
                min_delta = d;
                j0 = j_n[j];
                j_n_index = j;
            }
        }
        if (min_delta >= 0)
            return std::make_pair(x, j_basis);

        for (size_t i = 0; i < m; i++) col_to_replace[i] = A[i][j0];
        z_vec = matrix_vector_product<T>(B, col_to_replace);
        
        T min_omega = INF;
        size_t s;
        for(size_t j = 0; j < m;j++){
            if (z_vec[j] > 0){
                T om_i = x[j_basis[j]-1]/z_vec[j];
                if (om_i < min_omega){
                    min_omega = om_i;
                    s = j;
                }
            }
        }
        if (min_omega == INF){
            std::cout << "Unbounded\n";
            exit(0);
        }
        // new plan
        for (auto &&j : j_n){
            x[j] = 0;
        }
        x[j0] = min_omega;
        for(size_t i = 0; i < j_basis.size(); i++){
            x[j_basis[i]-1] -= min_omega*z_vec[i];
        }
        j_n[j_n_index] = j_basis[s] - 1;
        j_basis[s] = j0 + 1;
        c_basis[s] = c[j0];

        B = sherman_morrison<T>(B, col_to_replace, s);
    }
}

std::queue<size_t> intersection(const std::vector<size_t>& nums1, const std::vector<size_t>& nums2){
    std::unordered_set<size_t> set(nums2.cbegin(), nums2.cend());
    std::queue<size_t> intersect;
    for (const auto &n: nums1){
        if (set.erase(n) > 0)
            intersect.push(n);
    }
    return intersect;
}

template<typename T>
void transpose_matrix(std::vector<std::vector<T>>& A){
    for (size_t i = 0; i < A.size(); i++)
        for (size_t j = 0; j < A.size(); j++)
            std::swap(A[j][i], A[i][j]);
}

template<typename T>
std::vector<std::vector<T>> sherman_morrison2(const std::vector<std::vector<T>>& B, size_t k){
    std::vector<std::vector<T>> B_inv;
    std::vector<T> e(B.size());
    e[k] = 1;
    B_inv = sherman_morrison<T>(B, e, k);
    transpose_matrix<T>(B_inv);
    B_inv = sherman_morrison<T>(B_inv, e, k);
    B_inv.erase(B_inv.begin() + k);
    for(auto& row : B_inv) row.erase(row.begin() + k);
    return B_inv;
}

template<typename T>
std::pair<std::vector<T>, std::vector<size_t>> first_phase(std::vector<std::vector<T>>& A, std::vector<T> &b)
{
    size_t m = A.size(), n = A[0].size();
    std::vector<std::vector<T>> A_ext(m, std::vector<T>(m+n));
    std::vector<T> x(n+m), c_init(n);
    std::vector<size_t> j_basis0;

    for (size_t i = 0; i < m; i++){
        if (b[i] < 0){
            for(size_t j = 0; j<n; j++) A_ext[i][j] = -A[i][j];
            x[n+i] = -b[i];
        }
        else{
            x[n+i] = b[i];
            for(size_t j = 0; j<n; j++) A_ext[i][j] = A[i][j];
        }
        A_ext[i][n+i] = 1;
    }
    for (size_t i = 0; i < m; i++) c_init.push_back(-1);
    for (size_t i=n+1; i <= m+n; i++) j_basis0.push_back(i);
    
    std::pair<std::vector<T>, std::vector<size_t>> p = main_phase<T>(A_ext, c_init, x, j_basis0);
    for(size_t i = n; i < n+m; i++){
        if (abs(p.first[i]) > EPS){
            std::cout << "No solution" << '\n';
            exit(0);
        }
    }
    p.first.resize(n);
    //b,c
    std::queue<size_t> j_intersect = intersection(p.second, j_basis0);
    if(!j_intersect.empty()){
        std::vector<std::vector<T>> A_basis(m, std::vector<T>(m)), B;
        std::unordered_set<size_t> j_basis_set(p.second.begin(), p.second.end());
        std::vector<T> col_to_replace(m);

        for (size_t i = 0; i < m; i++)
            for(size_t j = 0; j < m; j++)
                A_basis[j][i] = A_ext[j][p.second[i]-1];
        B = invert_matrix<T>(A_basis);

        while(1)
        {
            auto it = std::find(p.second.begin(), p.second.end(), j_intersect.front());
            size_t k = it - p.second.begin();   
            size_t j0 = 0;
            for (size_t i = 1; i <= n; i++){
                if (j_basis_set.count(i) == 0){
                    T inner_pr = 0;
                    for (size_t j = 0; j < B[0].size(); j++){
                        inner_pr += B[k][j]*A_ext[j][i-1];
                    }
                    if (abs(inner_pr) > EPS){
                        j0 = i;
                        break;
                    }
                }
            }
            j_basis_set.erase(p.second[k]);
            j_intersect.pop();
            if (j0 != 0){
                p.second[k] = j0;
                if(j_intersect.empty()) break;
                j_basis_set.insert(j0);

                for (size_t i = 0; i < A_ext.size(); i++) col_to_replace[i] = A_ext[i][j0-1];
                B = sherman_morrison<T>(B, col_to_replace, k);
            }
            else{
                p.second.erase(it);
                A.erase(A.begin() + k);
                b.erase(b.begin() + k);
                if(j_intersect.empty()) break;

                B = sherman_morrison2<T>(B, k);
                A_ext.erase(A_ext.begin() + k);
                for(auto& row : A_ext) row.erase(row.begin() + *it-1);
                col_to_replace.resize(col_to_replace.size()-1);
            }
        }
    }
    return p;
}

int main() 
{ 
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    size_t m, n;
    std::cin >> m >> n;    
    std::vector<std::vector<double>> A(m, std::vector<double>(n));
    std::vector<double> c(n), b(m);
    for (size_t i = 0; i < m; i++){
        for (size_t j = 0; j < n; j++)
            std::cin >> A[i][j];
    }
    for (size_t i = 0; i < m; i++) std::cin >> b[i];
    for (size_t i = 0; i < n; i++) std::cin >> c[i];

    auto p0 = first_phase<double>(A, b);
    p0 = main_phase<double>(A, c, p0.first, p0.second);
    std::cout << "Bounded" << '\n';
    std::cout << std::setprecision(10) << std::fixed;
    for (auto &&el : p0.first) std::cout << el << ' ';
    std::cout<<'\n';
    return 0;
}
