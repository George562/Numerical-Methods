#include "tools.h"

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    ifstream input("input.txt");
    size_t n; input >> n;
    vector<double> X(n), foo(n);
    for (size_t i = 0; i < n; i ++) input >> X[i];
    for (size_t i = 0; i < n; i ++) input >> foo[i];
    double res = derivative(foo, X, 0.4, 1);
    cout << "res1 = " << res << '\n';
    res = derivative(foo, X, 0.4, 2);
    cout << "res2 = " << res << '\n';
    return 0;
}