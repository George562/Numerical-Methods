#include "tools.h"

double a        = 3;
double epsilon  = 0.0001;

int main() {
    makeGraph([](double x){ return 0.5 * sqrt(a * a - x * x); });
    makeGraph([](double x){ return -0.5 * sqrt(a * a - x * x); });
    makeGraph([](double x){ return (exp(x) + x) / a; });
    sf::Thread thread(ShowGraphics); thread.launch();

    Matrix<double> res, G(2, 2);
    cout << "l1 = "; cin >> G[0][0]; cout << "r1 = "; cin >> G[0][1];
    cout << "l2 = "; cin >> G[1][0]; cout << "r2 = "; cin >> G[1][1];

    cout << "Choose method: 1 - Newton, 2 - Iter\n> ";
    int choose; cin >> choose;
    MF foo(2, 1), J(2, 2);
    if (choose == 1) {
        foo[0][0] = [](Matrix<double> x) { return pow(x[0][0] / a, 2) + 4 * pow(x[1][0] / a, 2) - 1; };
        foo[1][0] = [](Matrix<double> x) { return a * x[1][0] - exp(x[0][0]) - x[0][0]; };

        J[0][0] = [](Matrix<double> x) { return 2 * x[0][0] / (a * a); };
        J[0][1] = [](Matrix<double> x) { return 8 * x[1][0] / (a * a); };
        J[1][0] = [](Matrix<double> x) { return - exp(x[0][0]) - 1; };
        J[1][1] = [](Matrix<double> x) { return a; };

        res = Newton(foo, J, G, epsilon);
        cout << res[0][0] << ' ' << res[1][0] << '\n';
    } else {
        foo[0][0] = [](Matrix<double> x) { return (a * x[1][0] - exp(x[0][0]) + 4 * x[0][0]) / 5; };
        foo[1][0] = [](Matrix<double> x) { return sqrt(a * a - x[0][0] * x[0][0]) / 2; };

        J[0][0] = [](Matrix<double> x) { return (4 - exp(x[0][0])) / 5; };
        J[0][1] = [](Matrix<double> x) { return a / 5; };
        J[1][0] = [](Matrix<double> x) { return -0.5 * x[0][0] / sqrt(a * a - x[0][0] * x[0][0]); };
        J[1][1] = [](Matrix<double> x) { return (double)0; };

        res = Iter(foo, J, G, epsilon);
        cout << res[0][0] << ' ' << res[1][0] << '\n';
    }

    thread.wait(); thread.terminate();
    return 0;
}
