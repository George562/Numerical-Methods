#include "../first/tools.h"
#include "../../lab1/TMA/TMA.hpp"
#include "../../Graphica.h"
using namespace std;

void CreateGraphic(vector<vector<double>> v) {
    vector<pair<double, double>> points;
    for (vector<double>& x: v) {
        addPoint(x[0], x[1]);
        points.emplace_back(x[0], x[1]);
    }
    makeGraphByPoints(points);
}

void print(size_t k, double heta_j, double y, double epsilon) {
    if (k == 0)
        cout << "-------------------------------------------------\n"
             << "|   j :   heta_j    :      y      : |F(heta_j)| |\n"
             << "------+-------------+-------------+--------------\n";
    cout << setprecision(4) << right << setfill(' ') << "| " << setw(3) << k << left;
    cout << setw(4) << " : " << setw(10) << heta_j;
    cout << setw(4) << " : " << setw(10) << y;
    cout << setw(4) << " : " << setw(10) << epsilon << " |\n";
    cout << "------+-------------+-------------+--------------\n";
}

vector<vector<double>> ShootingMethod(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double epsilon, double (*foo)(double)) {
    double heta_was = 10, heta_cur = 12, F_cur, F_was, heta_res;
    vector<vector<double>> res;

    res = Runge_Kutt_Method4(f, {x0[0], x0[1], heta_was}, a, b, h, foo);
    CreateGraphic(res);
    double y_was = res[res.size() - 1][1];
    F_was = res[res.size() - 1][1] - res[res.size() - 1][2] - 3.645;
    print(0, heta_was, y_was, abs(F_was));

    res = Runge_Kutt_Method4(f, {x0[0], x0[1], heta_cur}, a, b, h, foo);
    CreateGraphic(res);
    double y_cur = res[res.size() - 1][1];
    F_cur = res[res.size() - 1][1] - res[res.size() - 1][2] - 3.645;
    print(1, heta_cur, y_cur, abs(F_cur));

    for (size_t k = 2; abs(F_cur) > epsilon; k++) {
        heta_res = heta_cur - (heta_cur - heta_was) * F_cur / (F_cur - F_was);

        y_was = y_cur;
        heta_was = heta_cur;

        heta_cur = heta_res;
        res = Runge_Kutt_Method4(f, {x0[0], x0[1], heta_cur}, a, b, h, foo);
        CreateGraphic(res);
        y_cur = res[res.size() - 1][1];

        F_was = F_cur;
        F_cur = res[res.size() - 1][1] - res[res.size() - 1][2] - 3.645;
        print(k, heta_cur, y_cur, abs(F_cur));
    }

    return res;
}
vector<vector<double>> FiniteDifferenceMethod(double (*p)(double), double (*q)(double), double (*f)(double), double a, double b, double ya, double yb, double h) {
    vector<vector<double>> res;
    vector<double> qX, pX, fX;
    for (double x = a + h; x * sign(h) < b * sign(h); x += h) {
        pX.push_back(p(x));
        qX.push_back(q(x));
        fX.push_back(f(x));
    }
    size_t N = pX.size();

    Matrix<double> m(3, N), d(1, N);
    m[0][0] = 0; m[1][0] = 11; m[2][0] = 1;
    d[0][0] = 36.45;
    for (size_t i = 1; i < N - 1; i++) {
        m[0][i] = 1 - pX[i] * h / 2; m[1][i] = -2 + h * h * qX[i]; m[2][i] = 1 + pX[i] * h / 2;
        d[0][i] = h * h * fX[i];
    }
    m[0][N - 1] = 1 - pX[N - 1] * h / 2; m[1][N - 1] = -2 + h * h * qX[N - 1]; m[2][N - 1] = 0;
    d[0][N - 1] = h * h * fX[N - 1] - (1 + pX[N - 1] * h / 2) * yb;
    Matrix<double> y = TMAsolve(m, d);

    res.push_back({a, ya});
    for (size_t i = 0; i < N; i++) res.push_back({a + h * i, y[0][i]});
    res.push_back({b, yb});

    return res;
}
