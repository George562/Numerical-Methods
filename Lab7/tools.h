#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "../../Matrix.h"
#include "../Graphica.h"
using namespace std;

struct BorderCondition {
    double a, b;
    double (* func)(double, double);
    double operator() (double x, double y) { return func(x, y); }
};
using BC = BorderCondition;

int c;
void setC(double x) { c = x; }

namespace IterType {
    enum IterType : unsigned short int {
        simple,
        zaydel,
        relax
    };
};

Matrix<double> Iter(size_t Nx, size_t Ny, double Lx, double Ly, double epsilon,
                    BC leftCondition, BC rightCondition, BC downCondition, BC upCondition,
                    double (*f)(double, double), IterType::IterType flag, double omega = 1.) {
    double hx = Lx / Nx, hy = Ly / Ny;
    double hx2 = pow(hx, 2), hy2 = pow(hy, 2);

    Matrix<double> res(Nx + 1, Ny + 1, 0.35);

    double coef1 = - downCondition.a / hy;
    double coef2 =   downCondition.b + coef1;

    double coef3 = - leftCondition.a / hx;
    double coef4 =   leftCondition.b + coef3;

    double coef5 = upCondition.a / hy;
    double coef6 = upCondition.b + coef5;

    double coef7 = rightCondition.a / hx;
    double coef8 = rightCondition.b + coef7;

    for (size_t i =  1, j =  0; i <  Nx; i++) res[i][j] = (downCondition (i * hx, j * hy) + res[i]    [j + 1] * coef1) / coef2;
    for (size_t i =  0, j =  0; j <= Ny; j++) res[i][j] = (leftCondition (i * hx, j * hy) + res[i + 1][j]     * coef3) / coef4;
    for (size_t i =  1, j = Ny; i <  Nx; i++) res[i][j] = (upCondition   (i * hx, j * hy) + res[i]    [j - 1] * coef5) / coef6;
    for (size_t i = Nx, j =  0; j <= Ny; j++) res[i][j] = (rightCondition(i * hx, j * hy) + res[i - 1][j]     * coef7) / coef8;

    Matrix<double> was;
    double coef9 = 2 / hx2 + 2 / hy2;
    size_t counter = 0;
    do {
        counter++;
        was = res;
        for (size_t i = 1; i < Nx; i++) {
            for (size_t j = 1; j < Ny; j++) {
                if (flag == IterType::simple)
                    res[i][j] = ((was[i + 1][j] + was[i - 1][j]) / hx2 + (was[i][j + 1] + was[i][j - 1]) / hy2 + c * was[i][j] - hx2 * hy2 * f(i * hx, j * hy)) / coef9;
                else if (flag == IterType::zaydel)
                    res[i][j] = ((was[i + 1][j] + res[i - 1][j]) / hx2 + (was[i][j + 1] + res[i][j - 1]) / hy2 + c * was[i][j] - hx2 * hy2 * f(i * hx, j * hy)) / coef9;
                else if (flag == IterType::relax)
                    res[i][j] = omega * ((was[i + 1][j] + was[i - 1][j]) / hx2 + (was[i][j + 1] + was[i][j - 1]) / hy2 + c * was[i][j] - hx2 * hy2 * f(i * hx, j * hy)) / coef9
                                + (1 - omega) * was[i][j];
            }
        }
        if (flag == IterType::simple) {
            for (size_t i =  1, j =  0; i <  Nx; i++) res[i][j] = (downCondition (i * hx, j * hy) + was[i]    [j + 1] * coef1) / coef2;
            for (size_t i =  0, j =  0; j <= Ny; j++) res[i][j] = (leftCondition (i * hx, j * hy) + was[i + 1][j]     * coef3) / coef4;
            for (size_t i =  1, j = Ny; i <  Nx; i++) res[i][j] = (upCondition   (i * hx, j * hy) + was[i]    [j - 1] * coef5) / coef6;
            for (size_t i = Nx, j =  0; j <= Ny; j++) res[i][j] = (rightCondition(i * hx, j * hy) + was[i - 1][j]     * coef7) / coef8;
        } else {
            for (size_t i =  1, j =  0; i <  Nx; i++) res[i][j] = (downCondition (i * hx, j * hy) + res[i]    [j + 1] * coef1) / coef2;
            for (size_t i =  0, j =  0; j <= Ny; j++) res[i][j] = (leftCondition (i * hx, j * hy) + res[i + 1][j]     * coef3) / coef4;
            for (size_t i =  1, j = Ny; i <  Nx; i++) res[i][j] = (upCondition   (i * hx, j * hy) + res[i]    [j - 1] * coef5) / coef6;
            for (size_t i = Nx, j =  0; j <= Ny; j++) res[i][j] = (rightCondition(i * hx, j * hy) + res[i - 1][j]     * coef7) / coef8;
        }
    } while ((was - res).normL() > epsilon);
    cout << "counter = " << counter << '\n'; cout.flush();
    return res;
}
