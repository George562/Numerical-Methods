#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

size_t LDim;
vector<double> Lcoef, LXcopy;
size_t Lleft, Lright;

double (*Lagrange(double (*foo)(double), vector<double> X, double x, int dim))(double) {
    dim++;
    LDim = dim;
    Lcoef.resize(LDim);
    LXcopy = X;
    Lleft = 0; Lright = X.size();
    for (int i = 0; i < X.size() - dim; i++)
        if ((foo(X[Lleft]) + foo(X[Lright - 1])) / 2 > x) Lright--;
        else                                              Lleft++;
    for (size_t i = Lleft; i < Lright; i++) {
        double w = 1;
        for (size_t j = Lleft; j < Lright; j++)
            if (i != j)
                w *= X[i] - X[j];
        Lcoef[i] = foo(X[i]) / w;
    }
    return [](double x) {
        double ans = 0;
        for (size_t i = Lleft; i < Lright; i++) {
            double temp = Lcoef[i];
            for (size_t j = Lleft; j < Lright; j++)
                if (i != j)
                    temp *= x - LXcopy[j];
            ans += temp;
        }
        return ans;
    };
}

size_t NDim;
vector<double> Ncoef, NXcopy;
size_t Nleft, Nright;

double f(double (*foo)(double), vector<double> x, size_t left, size_t right) {
    if (right - left == 0) return foo(x[left]);
    return (f(foo, x, left, right - 1) - f(foo, x, left + 1, right)) / (x[left] - x[right]);
}

double (*Newton(double (*foo)(double), vector<double> X, double x, int dim))(double) {
    dim++;
    NDim = dim;
    Ncoef.resize(NDim);
    NXcopy = X;
    for (size_t i = 0; i < NDim; i++) Ncoef[i] = f(foo, X, 0, i);
    Nleft = 0; Nright = X.size();
    for (int i = 0; i < X.size() - dim; i++)
        if ((foo(X[Nleft]) + foo(X[Nright - 1])) / 2 > x) Nright--;
        else                                              Nleft++;
    return [](double x) {
        double ans = 0;
        for (size_t i = Nleft; i < Nright; i++) {
            double temp = Ncoef[i];
            for (size_t j = Nleft; j < i; j++)
                temp *= x - NXcopy[j];
            ans += temp;
        }
        return ans;
    };
}
