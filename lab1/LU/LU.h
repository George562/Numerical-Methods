#include "../../Matrix.h"
#include <cstdint>

template <typename T>
pair<Matrix<T>, vector<pair<size_t, size_t>>> LUdecompos(Matrix<T>& m) {
    if (m.dim.first != m.dim.second) {
        cerr << "LU decomposition working only for square matrix\n"; exit(1);
    }
    Matrix<T> LUm(m);
    vector<pair<size_t, size_t>> P;

    for (size_t j = 0, i, k; j < m.dim.first - 1; j++) {
        k = LUm.maxC(j, j);
        if (k != j) {
            LUm.swapL(j, k);
            P.emplace_back(j, k);
        }
        for (i = j + 1; i < m.dim.first; i++) {
            LUm[i][j] /= LUm[j][j];
            LUm.subL(j, i, LUm[i][j], j + 1);
        }
    }
    return {LUm, P};
}

template <typename T>
Matrix<T> LUsolve(Matrix<T>& m, Matrix<T>& b) {
    if (m.dim.first != m.dim.second || m.dim.first != b.dim.second) {
        cerr << "LU decomposition working only for square matrix\n"; exit(1);
    }
    Matrix<T> res(b);

    pair<Matrix<T>, vector<pair<size_t, size_t>>> resLU = LUdecompos(m);

    for (size_t n = 0, i, j; n < b.dim.first; n++) {
        for (pair<size_t, size_t>& x: resLU.second)
            swap(res[n][x.first], res[n][x.second]);
        for (i = 0; i < b.dim.second; i++) {
            for (j = 0; j < i; j++)
                res[n][i] -= res[n][j] * resLU.first[i][j];
        }
        cout.flush();
        for (i = b.dim.second - 1; i != UINT64_MAX; i--) {
            for (j = i + 1; j < b.dim.second; j++)
                res[n][i] -= res[n][j] * resLU.first[i][j];
            res[n][i] /= resLU.first[i][i];
        }
    }
    return res;
}
/* LUx=b
Lz=b
    b_1=z_1 -> z_1=b_1
    b_2=z_1*L_21 + z_2 => z_2=b_2-z_1*L_21
    b_k=sum_1^k(z_i*L_ki) => z_k=b_k-sum_1^{k-1}(z_i*L_ki)
Ux=z
    z_n=U_nn*x_n => x_n=z_n/U_nn
    z_{n-1}=U_{n-1}n*x_n+U_{n-1}{n-1}*x_{n-1} => x_{n-1}=(z_{n-1}-U_{n-1}n*x_n)/U_{n-1}{n-1}
    z_k=sum_k^n(U_ik) => x_k=(z_k-sum_{k+1}^n(U_ki*x_i))/U_kk
*/

template <typename T>
double detM(Matrix<T> m) {
    pair<Matrix<T>, vector<pair<size_t, size_t>>> resLU = LUdecompos(m);
    double res = 1;
    for (size_t i = 0; i < m.dim.first - 1; i++)
        res *= resLU.first[i][i];
    return res * ((resLU.second.size() % 2) ? -1 : 1);
}

template <typename T>
Matrix<T> inverseM(Matrix<T> m) {
    Matrix<T> e = E<T>(m.dim.second);
    return LUsolve(m, e);
}
