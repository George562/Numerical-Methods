#include <iostream>
#include <vector>
#include <fstream>
using namespace std;
#define diff(n) (foo[n + 1] - foo[n]) / (X[n + 1] - X[n])

double derivative(vector<double> foo, vector<double> X, double x, size_t p) {
    size_t Dim = X.size();
    if (x < X[0] || X[Dim - 1] < x) {
        cerr << "bruh\n"; exit(1);
    }
    size_t i;
    for (i = 0; X[i] < x; i++) {}
    i--;
    if (p == 1)
        return diff(i) + (2 * x - X[i] - X[i + 1]) * (diff(i + 1) - diff(i)) / (X[i + 2] - X[i]);
    else
        return 2 * (diff(i + 1) - diff(i)) / (X[i + 2] - X[i]);
}
