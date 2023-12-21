#include "tools.h"
#include "../../Graphica.h"

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    ifstream input("input.txt");
    size_t n; input >> n;
    vector<double> X(n), foo(n);
    for (size_t i = 0; i < n; i ++) input >> X[i];
    for (size_t i = 0; i < n; i ++) input >> foo[i];
    for (size_t i = 0; i < n; i ++) addPoint(X[i], foo[i]);
    double x; input >> x;
    double (*res)(double) = Spline(foo, X);
    cout << res(x) << '\n';

    SetDiapazon(1.05, 4.55);
    makeGraph(res);
    ShowGraphics();
    return 0;
}