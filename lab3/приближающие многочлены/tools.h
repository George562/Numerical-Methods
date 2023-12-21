#include "../../lab1/LU/LU.h"
#include "../../lab1/Iter/Iter.h"
size_t Dim;
vector<double> coef;

double (*Approximation(vector<double> foo, vector<double> X, size_t dim))(double) {
    Dim = dim + 1;
    coef.resize(Dim);
    Matrix<double> m(Dim, Dim), d(1, Dim);
    for (size_t i = 0; i < Dim; i++) {
        for (size_t j = 0; j < Dim; j++)
            if (i > j) m[j][i] = m[i][j];
            else
                for (size_t k = 0; k < X.size(); k++)
                    m[j][i] += pow(X[k], i + j);
        for (size_t j = 0; j < X.size(); j++)
            d[0][i] += foo[j] * pow(X[j], i);
    }
    Matrix<double> resM = LUsolve(m, d);
    for (size_t i = 0; i < Dim; i++)
        coef[i] = resM[0][i];
    return [](double x) {
        double res = coef[0];
        for (size_t i = 1; i < Dim; i++)
            res += coef[i] * pow(x, i);
        return res;
    };
}
