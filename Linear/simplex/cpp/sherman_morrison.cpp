#include <iostream> 
#include <iomanip>
#include <vector>


template<typename T>
std::vector<T> matrix_vector_product(const std::vector<std::vector<T>>& m, const std::vector<T>& v){
    std::vector<T> product_vec(v.size());
    for (int i = 0; i < v.size(); i++){
        for (int j = 0; j < v.size(); j++)
            product_vec[i] += v[j]*m[i][j];
    }
    return product_vec;
}

template<typename T>
std::vector<std::vector<T>> sherman_morrison(const std::vector<std::vector<T>>& B, const std::vector<T>& x, int k){
    std::vector<T> l_vec = matrix_vector_product<T>(B, x);
    T l_k = l_vec[k];
    if (l_k == 0){
        std::cout << "NO";
        exit(0);
    }
    l_vec[k] = -1;
    for (int i = 0; i < l_vec.size(); i++)
        l_vec[i] /= -l_k;
    std::vector<std::vector<T>> B_inv(B.size(), std::vector<T>(B.size()));
    for (int i = 0; i < B.size(); i++)
        if (i!=k){
            for (int j = 0; j < B.size(); j++)
                B_inv[i][j] = B[i][j] + l_vec[i]*B[k][j];
        }
        else{
            for (int j = 0; j < B.size(); j++)
                B_inv[i][j] = l_vec[k]*B[k][j];
        }
    return B_inv;
}

int main() 
{ 
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    int n, k;
    std::cin >> n >> k;
    --k;
    std::vector<std::vector<double>> A(n, std::vector<double>(n)), B(n, std::vector<double>(n));
    std::vector<double> x(n);

    for (int i = 0; i < n; i++) 
        for (int j = 0; j < n; j++) 
            std::cin >> A[i][j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) 
            std::cin >> B[i][j];
    for (int i = 0; i < n; i++)
        std::cin >> x[i];

    std::vector<std::vector<double>> B_inv = sherman_morrison<double>(B, x, k);
    std::cout << "YES" << '\n';
    std::cout << std::setprecision(6) << std::fixed;
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) 
            std::cout << B_inv[i][j] << ' ';
        std::cout << '\n';
    }
    return 0; 
}