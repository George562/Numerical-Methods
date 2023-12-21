#include "tools.h"

double foo(double x) { return pow(x / M_E, x) + pow(M_E / x, x); }

int main() {
    SetDiapazon(1.2, 4.5);
    double a = 1.2, b = 4.5, h = -0.1, epsilon = 10e-3;
    vector<vector<double>> res;
    char TheMethod;
    cout << "What TheMethod do you want to use?\n"
    "1 - ShootingMethod\n"
    "2 - FiniteDifferenceMethod\n";
    do { cout << ">>> "; cin >> TheMethod; } while (TheMethod < '1' || '2' < TheMethod);

    printing = false;
    if (TheMethod == '1') {
        vector<double> x = {b, 9.76690936686, 10};
        vector<double (*)(vector<double>)> f = {
            [](vector<double> x) { return x[2]; },
            [](vector<double> x) { return (x[2] + x[1] * x[0] * pow(log(x[0]), 3)) / (x[0] * log(x[0])); }
        };
        res = ShootingMethod(f, x, b, a, h, epsilon, foo);
    } else {
        res = FiniteDifferenceMethod([](double x) { return - 1. / (x * log(x)); },
        [](double x) { return - log(x) * log(x); }, [](double x) { return 0.; },
        1.2, 4.5, 3.04254890597, 9.76690936686, 0.1);
    }

    vector<pair<double, double>> points;
    for (vector<double>& x: res) {
        addPoint(x[0], x[1]);
        points.emplace_back(x[0], x[1]);
        // cout << "x = " << x[0] << "; y = " << x[1] << '\n';
    }
    makeGraphByPoints(points);
    makeGraph(foo, sf::Color::Red);
    ShowGraphics();
    return 0;
}