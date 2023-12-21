#include <iomanip>
#include "../../Matrix.h"
using namespace std;
#define sign(x) ((x > 0) ? 1 : -1)

bool printing = true;

Matrix<double> Butcher;

enum class Methods {
    Runge_Kutt_Method,
    ExplicitEulerMethod,
    MethodEuler_Cauchy,
    MethodAdams2,
    MethodAdams4
};

vector<vector<double>> Runge_Kutt_Method3(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double));
vector<vector<double>> Runge_Kutt_Method4(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double));
vector<vector<double>> ExplicitEulerMethod(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double));
vector<vector<double>> MethodEuler_Cauchy(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double));
vector<vector<double>> MethodAdams2(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double));
vector<vector<double>> MethodAdams4(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double));
vector<vector<double>> Method(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double), Methods);

void print(size_t k, vector<double> x, vector<double> Delta, double y, double epsilon) {
    if (k == 0)
        cout << "---------------------------------------------------------------------------------------------------------\n"
             << "|   k :      x      :      y      :      z      :      dy     :      dz     :     real    :   epsilon   |\n"
             << "------+-------------+-------------+-------------+-------------+-------------+-------------+--------------\n";
    cout << setprecision(4) << right << setfill(' ') << "| " << setw(3) << k << left;
    cout << setw(4) << " : " << setw(10) << x[0];
    cout << setw(4) << " : " << setw(10) << x[1];
    cout << setw(4) << " : " << setw(10) << x[2];
    cout << setw(4) << " : " << setw(10) << Delta[0];
    cout << setw(4) << " : " << setw(10) << Delta[1];
    cout << setw(4) << " : " << setw(10) << y;
    cout << setw(4) << " : " << setw(10) << epsilon << " |\n";
    cout << "------+-------------+-------------+-------------+-------------+-------------+-------------+--------------\n";
    // printf("| %3Ld | %.3lf | %.3lf | %.3lf | %.3lf | %.3lf | %.3lf | %.5lf |\n",
    //        k, x[0], x[1], x[2], Delta[0], Delta[1], y, epsilon);
}

vector<vector<double>> Method(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double), Methods method) {
    size_t Dim = f.size();
    vector<double> xk = x0, Delta(Dim);
    size_t k = 0;
    vector<vector<double>> res = { xk };
    if (method == Methods::MethodAdams4) {
        res = Runge_Kutt_Method4(f, x0, a, a + h * 3, h, foo);
        xk = res[3];
        k = 3;
    } else if (method == Methods::MethodAdams2) {
        res = Runge_Kutt_Method4(f, x0, a, a + h, h, foo);
        xk = res[1];
        k = 1;
    }
    for (; xk[0] * sign(h) < b * sign(h); k++) {
        switch (method) {
            case Methods::Runge_Kutt_Method: {
                size_t BDim = Butcher.dim.first;
                vector<vector<double>> K(Dim, vector<double>(BDim - 1, 0));
                for (size_t j = 0; j < BDim - 1; j++) {
                    vector<double> xj = xk;

                    xj[0] += Butcher[j][0] * h;
                    for (size_t i = 0; i < Dim; i++)
                        for (size_t k = 0; k < j; k++)
                            xj[i + 1] += Butcher[j][k + 1] * K[i][k];

                    for (size_t i = 0; i < Dim; i++) K[i][j] = h * f[i](xj);
                }
                for (size_t i = 0; i < Dim; i++) {
                    Delta[i] = 0;
                    for (size_t j = 0; j < BDim - 1; j++)
                        Delta[i] += Butcher[BDim - 1][j + 1] * K[i][j];
                }
                break;
            }

            case Methods::ExplicitEulerMethod:
                for (size_t i = 0; i < Dim; i++) Delta[i] = h * f[i](xk);
                break;

            case Methods::MethodEuler_Cauchy:
                for (size_t i = 0; i < Dim; i++) {
                    vector<double> xj = xk; xj[0] += h;
                    xj[i] = xk[i + 1] + h * f[i](xk);
                    Delta[i] = h * ( f[i](xk) + f[i](xj) ) / 2;
                }
                break;

            case Methods::MethodAdams2:
                for (size_t i = 0; i < Dim; i++)
                    Delta[i] = h * ( 3 * f[i](xk) - f[i](res[k - 1]) ) / 2;
                break;

            case Methods::MethodAdams4:
                for (size_t i = 0; i < Dim; i++)
                    Delta[i] = h * ( 55 * f[i](xk) - 59 * f[i](res[k - 1]) + 37 * f[i](res[k - 2]) - 9 * f[i](res[k - 3])) / 24;
                break;
        }
        if (printing) print(k, xk, Delta, foo(xk[0]), abs(foo(xk[0]) - xk[1]));

        xk[0] += h;
        for (size_t i = 0; i < Dim; i++) xk[i + 1] += Delta[i];
        res.push_back(xk);
    }
    if (printing) print(k, xk, Delta, foo(xk[0]), abs(foo(xk[0]) - xk[1]));
    return res;
}

vector<vector<double>> Runge_Kutt_Method2(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    Butcher.resize(3, 3);
    Butcher[0] = {    0,    0,    0 };
    Butcher[1] = {    1,    1,    0 };
    Butcher[2] = {    0,  0.5,  0.5 };
    return Method(f, x0, a, b, h, foo, Methods::Runge_Kutt_Method);
}
vector<vector<double>> Runge_Kutt_Method3(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    Butcher.resize(4, 4);
    Butcher[0] = {    0,    0,    0,    0 };
    Butcher[1] = { 1./3, 1./3,    0,    0 };
    Butcher[2] = { 2./3,    0, 2./3,    0 };
    Butcher[3] = {    0, 1./4,    0, 3./4 };
    return Method(f, x0, a, b, h, foo, Methods::Runge_Kutt_Method);
}
vector<vector<double>> Runge_Kutt_Method4(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    Butcher.resize(5, 5);
    Butcher[0] = {   0,    0,    0,    0,    0 };
    Butcher[1] = { 0.5,  0.5,    0,    0,    0 };
    Butcher[2] = { 0.5,    0,  0.5,    0,    0 };
    Butcher[3] = {   1,    0,    0,    1,    0 };
    Butcher[4] = {   0, 1./6, 1./3, 1./3, 1./6 };
    return Method(f, x0, a, b, h, foo, Methods::Runge_Kutt_Method);
}
vector<vector<double>> ExplicitEulerMethod(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    return Method(f, x0, a, b, h, foo, Methods::ExplicitEulerMethod);
}
vector<vector<double>> MethodEuler_Cauchy(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    return Method(f, x0, a, b, h, foo, Methods::MethodEuler_Cauchy);
}
vector<vector<double>> MethodAdams2(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    return Method(f, x0, a, b, h, foo, Methods::MethodAdams2);
}
vector<vector<double>> MethodAdams4(vector<double (*)(vector<double>)> f, vector<double> x0, double a, double b, double h, double (*foo)(double)) {
    return Method(f, x0, a, b, h, foo, Methods::MethodAdams4);
}