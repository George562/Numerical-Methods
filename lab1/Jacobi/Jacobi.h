#include "../../Matrix.h"

template <typename T>
T t(Matrix<T>& m) {
    T res = 0;
    for (size_t i = 0; i < m.dim.first; i++)
        for (size_t j = i + 1; j < m.dim.first; j++)
            res += m[i][j] * m[i][j];
    return pow(res, 0.5);
}

template <typename T>
pair<vector<T>, Matrix<T>> Jacobi(Matrix<T>& m, double epsilon) {
    size_t N = m.dim.second;
    for (size_t i = 0; i < N; i++)
        for (size_t j = i + 1; j < N; j++)
            if (m[i][j] != m[j][i]) {
                cerr << "Jacobi working only for simmetric matrix\n"; exit(1);
            }
    Matrix<T> alpha(m), u(N);
    Matrix<T> resM = E<T>(N);
    do {
        for (size_t i = 0; i < N; i++)
            for (size_t j = i + 1; j < N; j++) {
                u = E<T>(N);
                double phi = (alpha[i][i] == alpha[j][j]) ? M_PI / 4 : 0.5 * atan(2 * alpha[i][j] / (alpha[i][i] - alpha[j][j]));
                u[i][i] = cos(phi); u[i][j] = -sin(phi);
                u[j][i] = sin(phi); u[j][j] =  cos(phi);
                alpha = u.T() * alpha * u;
                resM = resM * u;
            }
    } while (t(alpha) > epsilon);
    vector<T> resV(N);
    for (size_t i = 0; i < N; i++) resV[i] = alpha[i][i];
    return {resV, resM};
}
