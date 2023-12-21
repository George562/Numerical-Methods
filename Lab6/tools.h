#include <iostream>
#include <cmath>
#include <vector>
#include "TMA.h"
#include <iomanip>
#include "../Graphica.h"
using namespace std;

double dxt0Func  (double x, double t) { return -exp(-x) * (cos(x) + sin(x)) ; }
double ddxt0Func (double x, double t) { return 2 * exp(-x) * sin(x)         ; }

double a, b, c;
void setABC(double x, double y, double z) { a = x; b = y; c = z; }

struct BorderCondition {
    double a, b;
    double (* func)(double, double);
    double operator() (double x, double y) { return func(x, y); }
};
using BC = BorderCondition;

Matrix<double> explicitMethod(size_t N, double L, double T,
                              BC leftCondition, BC rightCondition, BC t0Condition, BC dt0Condition,
                              double (*f)(double, double)) {
    double h = L / N;
    double h2 = pow(h, 2.);
    size_t K = T / (0.495 * h2);
    double tau = T / K;
    double tau2 = pow(tau, 2.);
    double sigma = tau2 / h2;

    double coef1 = a * a * sigma + b * sigma * h / 2;
    double coef2 = 2 + c * tau2 - 2 * a * a * sigma;
    double coef3 = a * a * sigma - b * sigma * h / 2;

    double coef6 = -leftCondition.a / h;
    double coef7 = rightCondition.a / h;

    double coef4 = leftCondition.b + coef6;
    double coef5 = rightCondition.b + coef7;

    Matrix<double> u(N + 1, K + 1);
    for (size_t j = 0; j <= N; j++) u[j][0] = t0Condition(j * h, 0.);
    for (size_t j = 0; j <= N; j++)
        u[j][1] = t0Condition(j * h, 0.)
                + tau * dt0Condition(j * h, 0.)
                + (tau2 / 2) * (a * a * ddxt0Func(j * h, 0.)
                              + b * dxt0Func(j * h, 0.)
                              + c * t0Condition(j * h, 0.)
                              + f(j * h, tau))
                ;

    for (size_t k = 1; k < K; k++) {
        for (size_t j = 1; j < N; j++)
            u[j][k + 1] = u[j + 1][k] * coef1
                        + u[j][k]     * coef2
                        + u[j - 1][k] * coef3
                        - u[j][k - 1]
                        + tau2 * f(j * h, k * tau);
        u[0][k + 1] = (leftCondition(0, (k + 1) * tau) + u[1][k + 1] * coef6) / coef4;
        u[N][k + 1] = (rightCondition(0, (k + 1) * tau) + u[N - 1][k + 1] * coef7) / coef5;
    }
    return u;
}

Matrix<double> implicitMethod(size_t N, double L, double T,
                              BC leftCondition, BC rightCondition, BC t0Condition, BC dt0Condition,
                              double (*f)(double, double)) {
    double h = L / N;
    double h2 = pow(h, 2.);
    size_t K = T / (0.495 * h2);
    double tau = T / K;
    double tau2 = pow(tau, 2.);
    double sigma = tau2 / h2;

    double coefA = -a * a * sigma + b * sigma * h / 2;
    double coefB = 1 + 2 * a * a * sigma - c * tau2;
    double coefC = -a * a * sigma - b * sigma * h / 2;

    Matrix<double> u(N + 1, K + 1);
    for (size_t j = 0; j <= N; j++) u[j][0] = t0Condition(j * h, 0.);
    for (size_t j = 0; j <= N; j++)
        u[j][1] = t0Condition(j * h, 0.)
                + tau * dt0Condition(j * h, 0.)
                + (tau2 / 2) * (a * a * ddxt0Func(j * h, 0.)
                              + b * dxt0Func(j * h, 0.)
                              + c * t0Condition(j * h, 0.)
                              + f(j * h, tau))
                ;

    Matrix<double> m(3, N + 1), d(1, N + 1);
    m[0][0] = 0;
    m[1][0] = leftCondition.b - leftCondition.a / h;
    m[2][0] = leftCondition.a / h;
    for (size_t j = 1; j < N; j++) {
        m[0][j] = coefA;
        m[1][j] = coefB;
        m[2][j] = coefC;
    }
    m[0][N] = -rightCondition.a / h;
    m[1][N] = rightCondition.b + rightCondition.a / h;
    m[2][N] = 0;
    for (size_t k = 1; k < K; k++) {
        for (size_t j = 1; j < N; j++)
            d[0][j] = 2 * u[j][k] - u[j][k - 1] + tau2 * f(j * h, (k + 1) * tau);
        d[0][0] = leftCondition(0, (k + 1) * tau);
        d[0][N] = rightCondition(0, (k + 1) * tau);
        u.copy(TMAsolve(m, d).T(), 0, k + 1);
    }
    return u;
}
