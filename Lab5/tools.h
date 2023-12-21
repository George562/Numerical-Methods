#include <iostream>
#include <cmath>
#include <vector>
#include "TMA.h"
#include <iomanip>
#include "../Graphica.h"
using namespace std;

struct BorderCondition {
    double a, b;
    double (* func)(double, double);
    double operator() (double x, double y) { return func(x, y); }
};
using BC = BorderCondition;

Matrix<double> explicitMethod(size_t N, double L, double T,
                              BC leftCondition, BC rightCondition, BC initialConditions,
                              double (*f)(double, double)) {
    double h = L / N;
    size_t K = T / (0.495 * pow(h, 2.));
    double tau = T / K;
    double sigma = tau / pow(h, 2.);

    Matrix<double> u(N + 1, K + 1);
    for (size_t j = 0; j <= N; j++) u[j][0] = initialConditions(j * h, 0.);

    for (size_t k = 0; k < K; k++) {
        for (size_t j = 1; j < N; j++)
            u[j][k + 1] = u[j][k] + sigma * (u[j + 1][k] - 2 * u[j][k] + u[j - 1][k]) + tau * f(j * h, k * tau);
        u[0][k + 1] = (leftCondition(0, (k + 1) * tau) - u[1][k + 1] * leftCondition.a / h) / (leftCondition.b - leftCondition.a / h);
        u[N][k + 1] = (rightCondition(0, (k + 1) * tau) + u[N - 1][k + 1] * rightCondition.a / h) / (rightCondition.b + rightCondition.a / h);
    }
    return u;
}

Matrix<double> implicitMethod(size_t N, double L, double T,
                              BC leftCondition, BC rightCondition, BC initialConditions,
                              double (*f)(double, double)) {
    double h = L / N;
    size_t K = T / (0.495 * pow(h, 2.));
    double tau = T / K;
    double sigma = tau / pow(h, 2.);

    Matrix<double> u(N + 1, K + 1);
    for (size_t j = 0; j <= N; j++) u[j][0] = initialConditions(j * h, 0.);

    Matrix<double> m(3, N + 1), d(1, N + 1);
    m[0][0] = 0;
    m[1][0] = leftCondition.b - leftCondition.a / h;
    m[2][0] = leftCondition.a / h;
    for (size_t j = 1; j < N; j++) {
        m[0][j] = sigma;
        m[1][j] = -(1 + 2 * sigma);
        m[2][j] = sigma;
    }
    m[0][N] = -rightCondition.a / h;
    m[1][N] = rightCondition.b + rightCondition.a / h;
    m[2][N] = 0;
    for (size_t k = 0; k < K; k++) {
        for (size_t j = 1; j < N; j++)
            d[0][j] = -u[j][k] - tau * f(j * h, (k + 1) * tau);
        d[0][0] = leftCondition(0, (k + 1) * tau);
        d[0][N] = rightCondition(0, (k + 1) * tau);
        u.copy(TMAsolve(m, d).T(), 0, k + 1);
    }
    return u;
}

Matrix<double> finiteDifferenceMethod(size_t N, double L, double T,
                                      BC leftCondition, BC rightCondition, BC initialConditions,
                                      double (*f)(double, double),
                                      double theta,
                                      int approximationType = 0) {
    double h = L / N;
    size_t K = T / (0.495 * pow(h, 2.));
    double tau = T / K;
    double sigma = tau / pow(h, 2.);

    Matrix<double> u(N + 1, K + 1);
    for (size_t j = 0; j <= N; j++) u[j][0] = initialConditions(j * h, 0.);

    Matrix<double> m(3, N + 1), d(1, N + 1);
    m[0][0] = 0;
    m[1][0] = leftCondition.b - leftCondition.a / h;
    if (approximationType == 0 || theta == 0) {
        m[2][0] = leftCondition.a / h;
    } else if (approximationType == 1) {
        m[2][0] = leftCondition.a / h - leftCondition.a * h / (2 * tau * theta);
    }
    for (size_t j = 1; j < N; j++) {
        m[0][j] = sigma * theta;
        m[1][j] = -(1 + 2 * sigma * theta);
        m[2][j] = sigma * theta;
    }
    if (approximationType == 0 || theta == 0) {
        m[0][N] = -rightCondition.a / h;
    } else if (approximationType == 1) {
        m[0][N] = -rightCondition.a / h + rightCondition.a * h / (2 * tau * theta);
    }
    m[1][N] = rightCondition.b + rightCondition.a / h;
    m[2][N] = 0;
    for (size_t k = 0; k < K; k++) {
        for (size_t j = 1; j < N; j++)
            d[0][j] = -(
                u[j][k] +
                tau * theta * f(j * h, (k + 1) * tau) +
                (1 - theta) * (sigma * (u[j + 1][k] - 2 * u[j][k] + u[j - 1][k]) + tau * f(j * h, k * tau))
            );
        if (approximationType == 0 || theta == 0) {
            d[0][0] = leftCondition(0, (k + 1) * tau);
            d[0][N] = rightCondition(0, (k + 1) * tau);
        } else if (approximationType == 1) {
            d[0][0] = leftCondition(0, (k + 1) * tau) + d[0][1] * leftCondition.a * h / (2 * tau * theta);
            d[0][N] = rightCondition(0, (k + 1) * tau) - d[0][N - 1] * rightCondition.a * h / (2 * tau * theta);
        }
        u.copy(TMAsolve(m, d).T(), 0, k + 1);
    }
    return u;
}

void print(Matrix<double> A) {
    for (int j = 0; j < A.dim.second; j++) cout << "+-----------";
    cout << "+\n";
    for (int i = 0; i < A.dim.first; i++) {
        cout << "| ";
        for (int j = 0; j < A.dim.second; j++) cout << setw(9) << to_string(A[i][j]) << ((j != A.dim.second - 1) ? " : " : "");
        cout << " |\n";
        for (int j = 0; j < A.dim.second; j++) cout << "+-----------";
        cout << "+\n";
    }
    cout << '\n';
}