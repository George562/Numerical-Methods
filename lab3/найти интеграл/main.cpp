#include "tools.h"

double foo(double x) { return x / (pow(x, 4) + 81); }

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);

    double X0 = 0, Xk = 2, h1 = 0.5, h2 = 0.25, ans = 0.0232346849766;
    double rect1 = rectangle(foo, X0, Xk, h1),  rect2 = rectangle(foo, X0, Xk, h2),
           trap1 = trapezoid(foo, X0, Xk, h1),  trap2 = trapezoid(foo, X0, Xk, h2),
           simp1 = Simpson(foo, X0, Xk, h1),    simp2 = Simpson(foo, X0, Xk, h2);
    cout << "rect1    = " << rect1 << '\n'
         << "error1   = " << abs(ans - rect1) << '\n' << '\n'
         << "rect2    = " << rect2 << '\n'
         << "error2   = " << abs(ans - rect2) << '\n' << '\n'
         << "trap1    = " << trap1 << '\n'
         << "error1   = " << abs(ans - trap1) << '\n' << '\n'
         << "trap2    = " << trap2 << '\n'
         << "error2   = " << abs(ans - trap2) << '\n' << '\n'
         << "Simpson1 = " << simp1   << '\n'
         << "error1   = " << abs(ans - simp1) << '\n' << '\n'
         << "Simpson2 = " << simp2   << '\n'
         << "error2   = " << abs(ans - simp2) << '\n';

    return 0;
}