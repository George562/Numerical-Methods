#include "QR.h"

double epsilon = 0.01;

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    Matrix<double> matrix((char*)"input.txt");
    auto [R, C] = QRsolve(matrix, epsilon);
    cout << "R = \n";
    for (size_t i = 0; i < R.size(); i++) cout << R[i] << '\n';
    cout << "C = \n";
    for (size_t i = 0; i < C.size(); i++)
        cout << C[i].first << " + " << C[i].second << "i\n";
    return 0;
}