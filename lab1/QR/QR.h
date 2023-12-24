#include "../../Matrix.h"

template <typename T>
T dist(pair<T, T>& left, pair<T, T>& right) {
    return pow(pow(left.first - right.first, 2) + pow(left.second - right.second, 2), 0.5);
}

template <typename T>
pair<Matrix<T>, Matrix<T>> QRdecompose(Matrix<T>& m) {
    size_t N = m.dim.second;
    Matrix<T> Q = E<T>(N), R(m), v(N, 1), H(N);

    for (size_t i = 0; i < N - 1; i++) {
        T temp = 0;
        for (size_t j = i; j < N; j++) {
            v[j][0] = R[j][i];
            temp = temp + R[j][i] * R[j][i];
        }
        v[i][0] = v[i][0] + ((R[i][i] < 0) ? -pow(temp, 0.5) : pow(temp, 0.5));
        H = E<T>(N) - (v * v.T()) * (T)(2.d / (v.T() * v)[0][0]);
        R = H * R;
        Q = Q * H;
        v[i][0] = 0;
    }
    return {Q, R};
}

template <typename T>
pair<vector<T>, vector<pair<T, T>>> QRsolve(Matrix<T>& m, double epsilon) {
    size_t N = m.dim.second;
    Matrix<T> A(m);
    vector<pair<T, T>> resC(N);
    vector<T> resR(N);
    bool flag;
    vector<bool> complex(N, false);
    do {
        flag = true;
        pair<Matrix<T>, Matrix<T>> qr = QRdecompose(A);
        A = qr.second * qr.first;
        for (size_t i = 0; i < N; i++) {
            T sum = (T)0;
            for (size_t j = i + 1; j < N; j++) sum = sum + A[j][i] * A[j][i];
            if (pow(sum, 0.5) > epsilon) {
                complex[i] = complex[i + 1] = true;
                T b = A[i][i] + A[i + 1][i + 1], c = A[i + 1][i] * A[i][i + 1];
                T d = pow(b, 2) - 4 * c;
                pair<T, T> sqr1 = {-b / 2.d,   pow(d, 0.5) / 2.d},
                           sqr2 = {-b / 2.d, - pow(d, 0.5) / 2.d};
                flag = flag && (max(dist(sqr1, resC[i]), dist(sqr2, resC[i + 1])) < epsilon);
                resC[i++] = sqr1;
                resC[i] = sqr2;
            } else {
                resR[i] = A[i][i];
                complex[i] = false;
            }
        }
    } while (!flag);
    size_t l = 0, k = 0;
    for (size_t i = 0; i < N; i++)
        if (complex[i]) resC[l++] = resC[i];
        else resR[k++] = resR[i];
    resR.resize(k);
    resC.resize(l);
    return {resR, resC};
}
