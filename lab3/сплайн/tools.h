#include <iostream>
#include <vector>
#include <fstream>
#include "../../lab1/TMA/TMA.hpp"
size_t Dim;
vector<double> Xcopy;
Matrix<double> coef;

double (*Spline(vector<double> foo, vector<double> X))(double) {
    Dim = X.size() - 1;
    coef.resize(Dim, Dim); coef[2][0] = 0;
    Xcopy = X;
    vector<double> h(Dim + 1, 0);
    for (size_t i = 1; i <= Dim; i++) h[i] = X[i] - X[i - 1];
    Matrix<double> m(3, Dim - 1), d(1, Dim - 1);
    for (size_t i = 2; i <= Dim; i++) {
        m[0][i - 2] = (i == 2) ? 0 : h[i - 1];
        m[1][i - 2] = 2 * (h[i] + h[i - 1]);
        m[2][i - 2] = (i == Dim) ? 0 : h[i];

        d[0][i - 2] = 3 * ((foo[i] - foo[i - 1]) / h[i] - (foo[i - 1] - foo[i - 2]) / h[i - 1]);
    }
    Matrix<double> c = TMAsolve(m, d);
    for (size_t i = 0; i < Dim; i++) {
        coef[0][i] = foo[i];
        if (i < Dim - 1) {
            coef[2][i + 1] = c[0][i];
            coef[1][i] = (foo[i + 1] - foo[i]) / h[i + 1] - h[i + 1] * (coef[2][i + 1] + 2 * coef[2][i]) / 3;
            coef[3][i] = (coef[2][i + 1] - coef[2][i]) / (3 * h[i + 1]);
        } else {
            coef[1][i] = (foo[i + 1] - foo[i]) / h[i + 1] - 2 * h[i + 1] * coef[2][i] / 3;
            coef[3][i] = -coef[2][i] / (3 * h[i + 1]);
        }
    }
    return [](double x) {
        if (x < Xcopy[0] || Xcopy[Dim] < x) {
            cerr << "bruh\n"; exit(1);
        }
        size_t i;
        for (i = 0; Xcopy[i] < x; i++) {}
        i--;
        return coef[0][i] + coef[1][i] * pow((x - Xcopy[i]), 1) +
               coef[2][i] * pow((x - Xcopy[i]), 2) + coef[3][i] * pow((x - Xcopy[i]), 3);
    };
}
