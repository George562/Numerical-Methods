#include <iomanip>
#include "TMA.h"
using namespace std;

struct BorderCondition {
    double a, b;
    double (* func)(double, double, double);
    double operator() (double x, double y, double t) { return func(x, y, t); }
};
using BC = BorderCondition;

double a = 1.;
void setA(double x) { a = x; }

vector<vector<vector<double>>> Method(size_t Nx, size_t Ny, size_t Nt, double Lx, double Ly, double Lt,
                                      BC leftCondition, BC rightCondition, BC downCondition, BC upCondition, BC initCondition,
                                      double (*f)(double, double, double), int type) {
    double hx = Lx / Nx, hy = Ly / Ny, tau = Lt / Nt;
    double hx2 = pow(hx, 2), hy2 = pow(hy, 2), ht_2 = tau / 2.;
    double sigmaX = a * tau / hx2, sigmaY = a * tau / hy2;

    vector<vector<vector<double>>> res(Nt + 1, vector<vector<double>>(Nx + 1, vector<double>(Ny + 1)));
    for (size_t i = 0; i < Nx + 1; i++)
        for (size_t j = 0; j < Ny + 1; j++)
            res[0][i][j] = initCondition(hx * i, hy * j, 0);

    Matrix<double> mx(3, Nx + 1), dx(1, Nx + 1), my(3, Ny + 1), dy(1, Ny + 1), u(Nx + 1, Ny + 1), r;

    double coef1 = - downCondition.a / hy;
    double coef2 =   downCondition.b + coef1;

    double coef3 = - leftCondition.a / hx;
    double coef4 =   leftCondition.b + coef3;

    double coef5 = upCondition.a / hy;
    double coef6 = upCondition.b + coef5;

    double coef7 = rightCondition.a / hx;
    double coef8 = rightCondition.b + coef7;

    // k + 1/2
    mx[0][0] = 0;
    mx[1][0] = leftCondition.b - leftCondition.a / hx;
    mx[2][0] = leftCondition.a / hx;
    for (size_t i = 1; i < Nx; i++) {
        mx[0][i] = sigmaX;
        mx[1][i] = -(2 + 2 * sigmaX);
        mx[2][i] = sigmaX;
    }
    mx[0][Nx] = -rightCondition.a / hx;
    mx[1][Nx] = rightCondition.b + rightCondition.a / hx;
    mx[2][Nx] = 0;
    // k + 1
    my[0][0] = 0;
    my[1][0] = downCondition.b - downCondition.a / hy;
    my[2][0] = downCondition.a / hy;
    for (size_t j = 1; j < Ny; j++) {
        my[0][j] = sigmaY;
        my[1][j] = -(2 + 2 * sigmaY);
        my[2][j] = sigmaY;
    }
    my[0][Ny] = -upCondition.a / hy;
    my[1][Ny] = upCondition.b + upCondition.a / hy;
    my[2][Ny] = 0;

    for (size_t k = 0; k < Nt; k++) {
        u.copy(Matrix(res[k]));
        // k + 1/2
        for (size_t j = 1; j < Ny; j++) {
            for (size_t i = 1; i < Nx; i++)
                dx[0][i] = -2 * res[k][i][j]
                           -sigmaY * (res[k][i][j + 1] - 2 * res[k][i][j] + res[k][i][j - 1]) * (1 - type)
                           -tau * f(i * hx, j * hy, (k + 0.5) * tau);
            dx[0][0] = leftCondition(0., j * hy, (k + 0.5) * tau);
            dx[0][Nx] = rightCondition(Lx, j * hy, (k + 0.5) * tau);
            u.copy(TMAsolve(mx, dx).T(), 0, j);
        }
        for (size_t j = 0; j <= Ny; j++) {
            res[k + 1][0] [j] = u[0] [j];
            res[k + 1][Nx][j] = u[Nx][j];
        }
        // k + 1
        for (size_t i = 1; i < Nx; i++) {
            for (size_t j = 1; j < Ny; j++)
                dy[0][j] = -2 * u[i][j]
                           -sigmaX * (u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) * (1 - type)
                           -tau * f(i * hx, j * hy, k * tau);
            dy[0][0] = downCondition(i * hx, 0., (k + 0.5) * tau);
            dy[0][Ny] = upCondition(i * hx, Ly, (k + 0.5) * tau);

            r = TMAsolve(my, dy);
            for (size_t j = 0; j <= Ny; j++)
                res[k + 1][i][j] = r[0][j];
        }
        res[k + 1][0] [0]  = (leftCondition (0., 0., (k + 0.5) * tau) + res[k + 1][1]     [0]  * coef3) / coef4;
        res[k + 1][0] [Ny] = (leftCondition (0., Ly, (k + 0.5) * tau) + res[k + 1][1]     [Ny] * coef3) / coef4;
        res[k + 1][Nx][0]  = (rightCondition(Lx, 0., (k + 0.5) * tau) + res[k + 1][Nx - 1][0]  * coef7) / coef8;
        res[k + 1][Nx][Ny] = (rightCondition(Lx, Ly, (k + 0.5) * tau) + res[k + 1][Nx - 1][Ny] * coef7) / coef8;
    }
    return res;
}

vector<vector<vector<double>>> AlternatingDirectionMethod(size_t Nx, size_t Ny, size_t Nt, double Lx, double Ly, double Lt,
                                    BC leftCondition, BC rightCondition, BC downCondition, BC upCondition, BC initCondition,
                                    double (*f)(double, double, double)) {
    return Method(Nx, Ny, Nt, Lx, Ly, Lt, leftCondition, rightCondition, downCondition, upCondition, initCondition, f, 0);
}
vector<vector<vector<double>>> FractionalStepsMethod(size_t Nx, size_t Ny, size_t Nt, double Lx, double Ly, double Lt,
                                    BC leftCondition, BC rightCondition, BC downCondition, BC upCondition, BC initCondition,
                                    double (*f)(double, double, double)) {
    return Method(Nx, Ny, Nt, Lx, Ly, Lt, leftCondition, rightCondition, downCondition, upCondition, initCondition, f, 1);
}
