#pragma once
#include "../Matrix.h"
#include <cstdint>
#include <thread>

////////////////////////////////////////////////////////////
// Class 
////////////////////////////////////////////////////////////

template <typename T> class LinearSolver {
private:
    Matrix<T> Iter(Matrix<T>& m, Matrix<T> b, double epsilon, bool flag);
public:
    // LU
    pair<Matrix<T>, vector<pair<size_t, size_t>>> LUdecompos(Matrix<T>& m);
    Matrix<T> LUsolve(Matrix<T>& m, Matrix<T>& b);
    double detM(Matrix<T>& m);
    Matrix<T> inverseM(Matrix<T>& m);

    // TMA
    bool CheckCondition(Matrix<T>& m);
    Matrix<T> TMAsolve(Matrix<T>& m, Matrix<T>& b);

    // Iter
    Matrix<T> SimpleIter(Matrix<T>& m, Matrix<T>& b, double epsilon);
    Matrix<T> Zaydel(Matrix<T>& m, Matrix<T>& b, double epsilon);

    // Jacobi
    T t(Matrix<T>& m);
    pair<vector<T>, Matrix<T>> Jacobi(Matrix<T>& m, double epsilon);

    // QR
    pair<Matrix<T>, Matrix<T>> QRdecompose(Matrix<T>& m);
    pair<vector<T>, vector<pair<T, T>>> QRsolve(Matrix<T>& m, double epsilon);
};

////////////////////////////////////////////////////////////
// Realization
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// LU
template <typename T>
pair<Matrix<T>, vector<pair<size_t, size_t>>> LinearSolver<T>::LUdecompos(Matrix<T>& m) {
    if (!m.isSquare()) {
        cerr << "LU decomposition working only for square matrix\n"; exit(1);
    }

    Matrix<T> LUm(m);
    vector<pair<size_t, size_t>> P;

    for (size_t j = 0, i; j < m.dim.first - 1; j++) {
        size_t k = LUm.maxC(j, j);
        if (k != j) {
            LUm.swapL(j, k);
            P.emplace_back(j, k);
        }
        #pragma omp parallel for shared(j, m, LUm) private(i)
        for (i = j + 1; i < m.dim.first; i++) {
            LUm[i][j] /= LUm[j][j];
            LUm.subL(j, i, LUm[i][j], j + 1);
        }
    }
    return {LUm, P};
}

template <typename T>
Matrix<T> LinearSolver<T>::LUsolve(Matrix<T>& m, Matrix<T>& b) {
    if (!m.isSquare()) {
        cerr << "LU decomposition working only for square matrix\n"; exit(1);
    }
    if (m.dim.first != b.dim.second) {
        cerr << "The dimension of the left side does not coincide with the dimension of the right side\n"; exit(1);
    }

    Matrix<T> res(b);
    pair<Matrix<T>, vector<pair<size_t, size_t>>> resLU = LUdecompos(m);

    size_t k;
    #pragma omp parallel for shared(resLU, res, b) private(k)
    for (k = 0; k < b.dim.first; k++) {
        for (size_t i = 0; i < resLU.second.size(); i++)
            swap(res[k][resLU.second[i].first], res[k][resLU.second[i].second]);
        for (size_t i = 0; i < res.dim.second; i++)
            for (size_t j = 0; j < i; j++)
                res[k][i] -= res[k][j] * resLU.first[i][j];
        for (size_t i = res.dim.second - 1; i != UINT64_MAX; i--) {
            for (size_t j = i + 1; j < res.dim.second; j++)
                res[k][i] -= res[k][j] * resLU.first[i][j];
            res[k][i] /= resLU.first[i][i];
        }
    }
    return res;
}

template <typename T>
double LinearSolver<T>::detM(Matrix<T>& m) {
    pair<Matrix<T>, vector<pair<size_t, size_t>>> resLU = LUdecompos(m);
    T res = resLU.first[0][0];
    for (size_t i = 1; i < m.dim.first; i++)
        res *= resLU.first[i][i];
    return res * ((resLU.second.size() % 2) ? -1 : 1);
}

template <typename T>
Matrix<T> LinearSolver<T>::inverseM(Matrix<T>& m) {
    Matrix<T> e = E<T>(m.dim.second);
    return LUsolve(m, e);
}
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// TMA
template <typename T>
bool LinearSolver<T>::CheckCondition(Matrix<T>& m) {
    bool res = true;
    int counter = 0;
    for (size_t i = 1; i < m.dim.second - 1; i++) {
        counter += (abs(m[1][i]) > abs(m[0][i]) + abs(m[2][i])) ? 1 : -1;
        res = (res && abs(m[1][i]) >= abs(m[0][i]) + abs(m[2][i]));
    }
    return res && (counter >= 0);
}

template <typename T>
Matrix<T> LinearSolver<T>::TMAsolve(Matrix<T>& m, Matrix<T>& b) {
    if (!CheckCondition(m)) cout << "Warning! for this matrix, the method is unstable!\n";
    size_t N = m.dim.second;
    Matrix<T> res(b);

    vector<T> P(N);
    P[0] = -m[2][0] / m[1][0]; P[N - 1] = 0;
    for (size_t i = 1; i < N - 1; i++)
        P[i] = -m[2][i] / (m[1][i] + m[0][i] * P[i - 1]);

    size_t k;
    #pragma omp parallel for shared(m, b, P, res, N) private(k)
    for (k = 0; k < b.dim.first; k++) {
        vector<T> Q(N);
        Q[0] = b[k][0] / m[1][0];
        for (size_t i = 1; i < N - 1; i++)
            Q[i] = (b[k][i] - m[0][i] * Q[i - 1]) / (m[1][i] + m[0][i] * P[i - 1]);
        Q[N - 1] = (b[k][N - 1] - m[0][N - 1] * Q[N - 2]) / (m[1][N - 1] + m[0][N - 1] * P[N - 2]);

        res[k][N - 1] = Q[N - 1];
        for (size_t i = N - 2; i != UINT64_MAX; i--)
            res[k][i] = res[k][i + 1] * P[i] + Q[i];
    }
    return res;
}
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Iter
template <typename T>
Matrix<T> LinearSolver<T>::Iter(Matrix<T>& m, Matrix<T> b, double epsilon, bool flag) {
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
    #pragma omp parallel for shared(res, N, b, epsilon, flag, alpha, normAlpha) private(k, was)
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
Matrix<T> LinearSolver<T>::SimpleIter(Matrix<T>& m, Matrix<T>& b, double epsilon) {
    return Iter(m, b, epsilon, true);
}
template <typename T>
Matrix<T> LinearSolver<T>::Zaydel(Matrix<T>& m, Matrix<T>& b, double epsilon) {
    return Iter(m, b, epsilon, false);
}
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Jacobi
template <typename T>
T LinearSolver<T>::t(Matrix<T>& m) {
    T res = 0;
    for (size_t i = 0; i < m.dim.first; i++)
        for (size_t j = i + 1; j < m.dim.first; j++)
            res += m[i][j] * m[i][j];
    return pow(res, 0.5d);
}

template <typename T>
pair<vector<T>, Matrix<T>> LinearSolver<T>::Jacobi(Matrix<T>& m, double epsilon) {
    size_t N = m.dim.second;
    if (!m.isSymmetrical()) { cerr << "Jacobi working only for simmetric matrix\n"; exit(1); }

    Matrix<T> alpha(m), u = E<T>(N), uT = E<T>(N), resM = E<T>(N);
    double phi;
    do {
        for (size_t i = 0; i < N; i++)
            for (size_t j = i + 1; j < N; j++) {
                phi = (alpha[i][i] == alpha[j][j]) ? M_PI_4 : 0.5 * atan(2 * alpha[i][j] / (alpha[i][i] - alpha[j][j]));
                u[i][i] = cos(phi); u[i][j] = -sin(phi);
                u[j][i] = sin(phi); u[j][j] =  cos(phi);

                uT[i][i] =  cos(phi); uT[i][j] = sin(phi);
                uT[j][i] = -sin(phi); uT[j][j] = cos(phi);

                
                std::thread first([&](){ alpha = uT * alpha * u; }),
                            second([&](){ resM = resM * u; });
                first.join();
                second.join();

                u[i][i] = T(1); u[i][j] = T(0);
                u[j][i] = T(0); u[j][j] = T(1);
                uT[i][i] = T(1); uT[i][j] = T(0);
                uT[j][i] = T(0); uT[j][j] = T(1);
            }
    } while (t(alpha) > epsilon);
    vector<T> resV(N);
    for (size_t i = 0; i < N; i++) resV[i] = alpha[i][i];
    return {resV, resM};
}
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// QR
template <typename T>
pair<Matrix<T>, Matrix<T>> LinearSolver<T>::QRdecompose(Matrix<T>& m) {
    size_t N = m.dim.second;
    Matrix<T> Q = E<T>(N), R(m), v(N, 1), H(N);

    for (size_t i = 0; i < N - 1; i++) {
        T temp = 0;
        for (size_t j = i; j < N; j++) {
            v[j][0] = R[j][i];
            temp = temp + R[j][i] * R[j][i];
        }
        v[i][0] += ((R[i][i] < 0) ? -pow(temp, 0.5) : pow(temp, 0.5));
        H = E<T>(N) - (v * v.T()) * (T)(2.d / (v.T() * v)[0][0]);
        std::thread first([&](){ R = H * R; }),
                    second([&](){ Q = Q * H; });
        v[i][0] = 0;
        first.join();
        second.join();
    }
    return {Q, R};
}

template <typename T>
pair<vector<T>, vector<pair<T, T>>> LinearSolver<T>::QRsolve(Matrix<T>& m, double epsilon) {
    size_t N = m.dim.second;
    Matrix<T> A(m);
    vector<pair<T, T>> resC(N);
    vector<T> resR(N);
    bool flag;
    size_t indexC, indexR;
    do {
        indexC = 0; indexR = 0;
        flag = true;
        pair<Matrix<T>, Matrix<T>> qr = QRdecompose(A);
        A = qr.second * qr.first;
        for (size_t i = 0; i < N; i++) {
            T sum = (T)0;
            for (size_t j = i + 1; j < N; j++) sum = sum + A[j][i] * A[j][i];
            if (pow(sum, 0.5) > epsilon) {
                T b = A[i][i] + A[i + 1][i + 1], c = A[i + 1][i] * A[i][i + 1];
                T d = pow(b, 2) - 4 * c;
                pair<T, T> sqr1 = {-b / 2.d,   pow(d, 0.5) / 2.d},
                           sqr2 = {-b / 2.d, - pow(d, 0.5) / 2.d};
                flag = flag && (max(hypot(sqr1.first - resC[indexC].first, sqr1.second - resC[indexC].second),
                                    hypot(sqr2.first - resC[indexC + 1].first, sqr2.second - resC[indexC + 1].second)) < epsilon);
                resC[indexC++] = sqr1;
                resC[indexC++] = sqr2;
                i += 1;
            } else {
                resR[indexR++] = A[i][i];
            }
        }
    } while (!flag);
    resR.resize(indexR);
    resC.resize(indexC);
    return {resR, resC};
}
////////////////////////////////////////////////////////////