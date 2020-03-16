#include <iostream> 
#include <iomanip>
#include <vector>
#include <algorithm>

#define INF 1e16

template<typename T>
std::vector<T> matrix_vector_product(const std::vector<std::vector<T>>& m, const std::vector<T>& v){
    std::vector<T> product_vec(m.size());
    for (int i = 0; i < m.size(); i++){
        for (int j = 0; j < v.size(); j++)
            product_vec[i] += m[i][j]*v[j];
    }
    return product_vec;
}

template<typename T>
std::vector<T> vector_matrix_product(const std::vector<T>& v, const std::vector<std::vector<T>>& m){
    std::vector<T> product_vec(m[0].size());
    for (int i = 0; i < m[0].size(); i++){
        for (int j = 0; j < v.size(); j++)
            product_vec[i] += v[j]*m[j][i];
    }
    return product_vec;
}

template<typename T>
std::vector<std::vector<T>> sherman_morrison(const std::vector<std::vector<T>>& B, const std::vector<T>& x, int k){
    std::vector<T> l_vec = matrix_vector_product<T>(B, x);
    T l_k = l_vec[k];
    // if (l_k == 0)
    l_vec[k] = -1;
    for (auto &el : l_vec)
        el /= -l_k;
    std::vector<std::vector<T>> B_inv(B.size(), std::vector<T>(B.size()));
    for (int i = 0; i < B.size(); i++)
        if (i != k){
            for (int j = 0; j < B.size(); j++)
                B_inv[i][j] = B[i][j] + l_vec[i]*B[k][j];
        }
        else{
            for (int j = 0; j < B.size(); j++)
                B_inv[i][j] = l_vec[k]*B[k][j];
        }
    return B_inv;
}

template<typename T>
std::vector<std::vector<T>> invert_matrix(std::vector<std::vector<T>>& A){
    std::vector<std::vector<T>> E(A.size(), std::vector<T>(A.size()));
    for (int i=0; i < E.size(); i++) ++E[i][i];

    for (int i = 0; i < A.size(); i++){
        int max_row = i;
        for (int j = i+1; j < A.size(); j++)
            if (abs(A[j][i]) > abs(A[max_row][i]))
                max_row = j;
        if (max_row != i){
            A[i].swap(A[max_row]);
            E[i].swap(E[max_row]);
        }
        for(int j = i+1; j < A.size(); j++){
            T c = A[j][i]/A[i][i];
            for(int k = 0; k < i; k++)
                E[j][k] -= E[i][k]*c;
            for(int k = i; k < A.size(); k++){
                A[j][k] -= A[i][k]*c;
                E[j][k] -= E[i][k]*c;
            }
        }
    }
    for(int i = A.size()-1; i >= 0; i--){
        T c = A[i][i];
        for(int j = 0; j < i; j++)
            for(int k = E.size()-1; k >= 0; k--)
                E[j][k] -= E[i][k]*A[j][i]/c;
        for(int j = 0; j < E.size(); j++)
            E[i][j] /= c;
    }
    return E;
}

template<typename T>
std::pair<std::vector<T>, std::vector<int>> main_phase(const std::vector<std::vector<T>>& A,
    const std::vector<T>& b, const std::vector<T>& c, std::vector<T> x, std::vector<int> j_basis)
{
    int m = A.size(), n = A[0].size();
    std::vector<std::vector<T>> A_basis(m, std::vector<double>(m)), B;
    std::vector<double> c_basis(m);
    std::vector<int> j_n;

    for (int i = 0; i < m; i++)
        c_basis[i] = c[j_basis[i]-1];
    for (int i=0; i < m; i++)
        for(int j=0; j < m; j++)
            A_basis[j][i] = A[j][j_basis[i]-1];

    std::vector<int> j_basis_sorted;
    std::copy(j_basis.begin(), j_basis.end(), std::back_inserter(j_basis_sorted));
    std::sort(j_basis_sorted.begin(), j_basis_sorted.end());
    for (int i=0, j=0; i < n; i++){
        if(j_basis_sorted[j]-1 == i && j < m) ++j;
        else j_n.push_back(i);
    }
    B = invert_matrix<T>(A_basis);
    std::vector<T> u_vec, uA_vec, col_to_replace(m), z_vec;

    while (1)
    {
        u_vec = vector_matrix_product<T>(c_basis, B);
        uA_vec = vector_matrix_product<T>(u_vec, A);
        T min_delta = INF;
        int j0, j_n_index;
        for (int j=0; j < j_n.size(); j++){
            T d = uA_vec[j_n[j]] - c[j_n[j]];
            if (d < min_delta){
                min_delta = d;
                j0 = j_n[j];
                j_n_index = j;
            }
        }
        if (min_delta >= 0)
            return std::make_pair(x, j_basis);

        for (int i = 0; i < m; i++) col_to_replace[i] = A[i][j0];
        z_vec = matrix_vector_product(B, col_to_replace);
        
        T min_omega = INF;
        int s;
        for(int j=0; j < m;j++){
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
        for(int i=0; i < j_basis.size(); i++){
            x[j_basis[i]-1] -= min_omega*z_vec[i];
        }
        j_n[j_n_index] = j_basis[s] - 1;
        j_basis[s] = j0 + 1;
        c_basis[s] = c[j0];

        B = sherman_morrison<T>(B, col_to_replace, s);
    }
}

int main() 
{ 
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    int m, n;
    std::cin >> m >> n;
    
    std::vector<std::vector<double>> A(m, std::vector<double>(n));
    std::vector<double> x(n), c(n), b(m);
    std::vector<int> j_basis(m);
    
    for (int i = 0; i < m; i++) 
        for (int j = 0; j < n; j++) 
            std::cin >> A[i][j];
    for (int i = 0; i < m; i++) std::cin >> b[i];
    for (int i = 0; i < n; i++) std::cin >> c[i];
    for (int i = 0; i < n; i++) std::cin >> x[i];
    for (int i = 0; i < m; i++) std::cin >> j_basis[i];

    std::pair<std::vector<double>, std::vector<int>> p = main_phase<double>(A, b, c, x, j_basis);
    for (auto &&el : A) {for (auto &&el1 : el) std::cout << el1 << ' '; std::cout<<'\n';}
    std::cout << "Bounded\n";
    std::cout << std::setprecision(10) << std::fixed;
    for (auto &&el : p.first) std::cout << el << ' ';
    std::cout << '\n';
    return 0;
}