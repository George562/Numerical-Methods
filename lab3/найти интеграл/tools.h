#include <iostream>
#include <cmath>
using namespace std;

double rectangle(double (*foo)(double), double x0, double xk, double h) {
    size_t Dim = (xk - x0) / h;
    double res = 0;
    for (size_t i = 0; i <= Dim; i++)
        res += foo(x0 + h * double(i * 2 + 1) / 2);
    return res * h;
}

double trapezoid(double (*foo)(double), double x0, double xk, double h) {
    size_t Dim = (xk - x0) / h;
    double res = (foo(x0) + foo(xk)) / 2;
    for (size_t i = 1; i < Dim; i++)
        res += foo(x0 + h * i);
    return res * h;
}

double Simpson(double (*foo)(double), double x0, double xk, double h) {
    if ((xk - x0) / (2 * h) != int((xk - x0) / (2 * h))) {
        cerr << "wrong interval or legth of step\n"; exit(1);
    }
    size_t Dim = (xk - x0) / h;
    double res = foo(x0) - foo(xk);
    size_t i;
    for (i = 1; i <= Dim; i++)
        res += 4 * foo(x0 + h * (double(i) - 0.5)) + 2 * foo(x0 + h * i);
    return res * h / 6;
}
