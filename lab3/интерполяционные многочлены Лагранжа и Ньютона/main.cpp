#include "tools.h"
#include "../../Graphica.h"

double foo(double x) { return tan(M_PI_2 - x) + x; }

int main() {
    vector<double> X1 = {M_PI_2 / 8, M_PI_2 / 4, M_PI_2 * 3 / 8, M_PI_2 / 2},
                   X2 = {M_PI_2 / 8, M_PI_2 / 4, M_PI_2 * 3 / 8, M_PI_2 / 2};
    SetDiapazon(X1[0], X1[3]);
    makeGraph(foo);
    for (size_t i = 0; i < X1.size(); i++) addPoint(X1[i], foo(X1[i]));
    double x = M_PI_2 * 3 / 16;

    int n; cout << "n for Lagrange = "; cin >> n;
    double (*ans1)(double) = Lagrange(foo, X1, x, n);
    double val1 = ans1(x), real1 = foo(x);
    cout << "val1 = " << val1 << "\nAbsolute interpolation error is: " << abs(val1 - real1) << '\n';
    makeGraph(ans1, sf::Color::Green);

    cout << "n for Newton = "; cin >> n;
    double (*ans2)(double) = Newton(foo, X2, x, n);
    double val2 = ans2(x), real2 = foo(x);
    cout << "val2 = " << val2 << "\nAbsolute interpolation error is: " << abs(val2 - real2) << '\n';
    makeGraph(ans2, sf::Color::Red);

    cout << "Green - Lagrange\nRed - Newton\n";
    ShowGraphics();
    return 0;
}