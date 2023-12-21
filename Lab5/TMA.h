#include "../Matrix.h"
#include <cstdint>

template <typename T>
bool CheckCondition(Matrix<T>& m) {
    bool res = true;
    int counter = 0;
    for (size_t i = 1; i < m.dim.second - 1; i++) {
        counter += (abs(m[1][i]) > abs(m[0][i]) + abs(m[2][i])) ? 1 : -1;
        res = (res && // предыдущее значение
              abs(m[1][i]) >= abs(m[0][i]) + abs(m[2][i])); // проверка качель
    }
    return res && (counter >= 0);
}

template <typename T>
Matrix<T> TMAsolve(Matrix<T>& m, Matrix<T>& d) {
    if (!CheckCondition(m)) cout << "Warning! for this matrix, the method is unstable!\n";
    size_t N = m.dim.second;
    Matrix<T> res(d);

    vector<T> P(N), Q(N);
    P[0] = -m[2][0] / m[1][0]; P[N - 1] = 0;
    for (size_t i = 1; i < N - 1; i++)
        P[i] = -m[2][i] / (m[1][i] + m[0][i] * P[i - 1]);

    for (size_t n = 0; n < d.dim.first; n++) {
        Q[0] = d[n][0] / m[1][0];
        for (size_t i = 1; i < N - 1; i++)
            Q[i] = (d[n][i] - m[0][i] * Q[i - 1]) / (m[1][i] + m[0][i] * P[i - 1]);
        Q[N - 1] = (d[n][N - 1] - m[0][N - 1] * Q[N - 2]) / (m[1][N - 1] + m[0][N - 1] * P[N - 2]);

        res[n][N - 1] = Q[N - 1];
        for (size_t i = N - 2; i != UINT64_MAX; i--)
            res[n][i] = res[n][i + 1] * P[i] + Q[i];
    }
    return res;
}