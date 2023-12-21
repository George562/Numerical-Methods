#include "tools.h"
#include "../../Graphica.h"

double foo(double x) { return  1 + log(abs(x)); }

int main() {
    vector<vector<double>> res;
    vector<vector<vector<double>> (*)(vector<double (*)(vector<double>)>, vector<double>, double, double, double, double (*)(double))> funcs = {
        Runge_Kutt_Method2,
        Runge_Kutt_Method3,
        Runge_Kutt_Method4,
        ExplicitEulerMethod,
        MethodEuler_Cauchy,
        MethodAdams2,
        MethodAdams4
    };
    double a = 1, b = 2, h = 0.1;
    vector<double> x0 = {a, 1, 1};
    vector<double (*)(vector<double>)> f(2);
    f[0] = [](vector<double> x) { return x[2]; };
    f[1] = [](vector<double> x) { return - x[2] / x[0]; };

    char TheMethod;
    cout << "What TheMethod do you want to use?\n"
    "1 - Runge_Kutt_Method2\n"
    "2 - Runge_Kutt_Method3\n"
    "3 - Runge_Kutt_Method4\n"
    "4 - ExplicitEulerMethod\n"
    "5 - MethodEuler_Cauchy\n"
    "6 - MethodAdams2\n"
    "7 - MethodAdams4\n";
    do { cout << ">>> "; cin >> TheMethod; } while (TheMethod < '1' || '7' < TheMethod);

    char reply;
    do { cout << "Want to use the standard value h = 0.1?\ny/n: "; cin >> reply;
    } while (reply != 'n' && reply != 'y');
    if (reply == 'n') { cout << "Enter a new h value: "; cin >> h; }

    res = funcs[TheMethod - '1'](f, x0, a, b, h, foo);
    vector<pair<double, double>> points;
    for (vector<double>& x: res) {
        addPoint(x[0], x[1]);
        points.emplace_back(x[0], x[1]);
    }
    SetDiapazon(a, b);
    makeGraphByPoints(points);
    makeGraph(foo, sf::Color::Red);
    ShowGraphics();
    return 0;
}