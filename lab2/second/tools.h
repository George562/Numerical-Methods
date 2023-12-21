#include "../../Graphica.h"
#include "../../Matrix.h"
#include <iostream>

using MF = Matrix<double (*)(Matrix<double>)>;

template <typename T>
double trashDetM(Matrix<T> m) {
    if (m.dim.first == 1) return m[0][0];
    if (m.dim.first == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    double res = 1;
    return res;
}

Matrix<double> subsM(MF m, Matrix<double> x) {
    Matrix<double> res(m.dim);
    for (size_t i = 0; i < m.dim.first; i++)
        for (size_t j = 0; j < m.dim.second; j++)
            res[i][j] = m[i][j](x);
    return res;
}

Matrix<double> Newton(MF foo, MF J, Matrix<double> G, double epsilon) {
    Matrix<double> x(G.dim.first, 1);
    for (size_t i = 0; i < x.dim.first; i++)
        x[i][0] = (G[i][1] + G[i][0]) / 2;

    Matrix<double> was;
    size_t counter = 0;
    do {
        counter++;
        was = x;
        double detJ = trashDetM(subsM(J, x));
        for (size_t i = 0; i < x.dim.first; i++) {
            MF Ax(J);
            for (size_t j = 0; j < x.dim.first; j++) Ax[j][i] = foo[j][0];
            x[i][0] = x[i][0] - trashDetM(subsM(Ax, was)) / detJ;
        }
    } while ((x - was).normC() > epsilon);
    cout << "counter = " << counter << '\n';
    return x;
}

Matrix<double> Iter(MF foo, MF foo1, Matrix<double> G, double epsilon) {
    size_t N = foo.dim.first;
    double q = 0;
    Matrix<double> was, x(N, 1);
    for (size_t i = 0; i < N; i++) x[i][0] = G[i][0];
    while (x[N - 1][0] < G[N - 1][1] > 0) {
        for (size_t i = 0; i < N; i++) {
            double temp = 0;
            for (size_t j = 0; j < N; j++) temp += abs(foo1[i][j](x));
            q = max(q, temp);
        }
        x[0][0] += (G[0][1] - G[0][0]) / 100;
        for (size_t i = 0; i < N - 1; i++)
            if (x[i][0] > G[i][1]) {
                x[i][0] = G[i][0];
                x[i + 1][0] += (G[i + 1][1] - G[i + 1][0]) / 100;
            }
    }
    cout << "q = " << q << '\n';
    for (size_t i = 0; i < N; i++) x[i][0] = (G[i][1] + G[i][0]) / 2;
    size_t counter = 0;
    do {
        counter++;
        was = x;
        x = subsM(foo, x);
    } while ((x - was).normC() > epsilon);
    cout << "counter = " << counter << '\n';
    return x;
}
