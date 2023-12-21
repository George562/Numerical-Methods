#include "../../Matrix.h"

template <typename T>
Matrix<T> Iter(Matrix<T>& m, Matrix<T> b, double epsilon, bool flag) {
    size_t N = m.dim.second;
    Matrix<T> alpha(m);

    T coef;
    for (size_t i = 0; i < N; i++) {
        coef = alpha[i][i];
        alpha[i][i] = 0;
        for (size_t j = 0; j < N; j++) {
            alpha[i][j] /= -coef;
        }
        for (size_t n = 0; n < b.dim.first; n++) {
            b[n][i] /= coef;
        }
    }

    double normAlpha = alpha.normC();
    Matrix<T> res(b.dim);

    Matrix<T> was;
    size_t k;
    for (k = 0; k < b.dim.first; k++) {
        Matrix<T> cur(b[k], N, 1), betha(b[k], N, 1);
        do {
            was = cur;
            if (flag) {
                cur = betha + alpha * cur;
            } else {
                for (size_t i = 0; i < N; i++) {
                    cur[i][0] = betha[i][0];
                    for (size_t j = 0; j < i; j++) { 
                        cur[i][0] += alpha[i][j] * cur[j][0];
                    }
                    for (size_t j = i; j < N; j++) { 
                        cur[i][0] += alpha[i][j] * was[j][0];
                    }
                }
            }
        } while ((cur - was).normC() * (normAlpha / (1 - normAlpha)) > epsilon);
        for (size_t i = 0; i < N; i++) {
            res[k][i] = cur[i][0];
        }
    }
    return res;
}

template <typename T>
Matrix<T> SimpleIter(Matrix<T>& m, Matrix<T>& b, double epsilon) {
    return Iter(m, b, epsilon, true);
}
template <typename T>
Matrix<T> Zaydel(Matrix<T>& m, Matrix<T>& b, double epsilon) {
    return Iter(m, b, epsilon, false);
}