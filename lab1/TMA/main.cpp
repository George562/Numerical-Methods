#include "TMA.hpp"

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    Matrix<double> matrix((char*)"inputABC.txt");
    Matrix<double> b((char*)"inputD.txt");
    Matrix<double> answer = TMAsolve(matrix, b);
    answer.recordM((char*)"answer.txt");
    return 0;
}