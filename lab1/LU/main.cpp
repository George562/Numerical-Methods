#include "LU.h"
#include <ctime>

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    Matrix<double> matrix((char*)"inputM.txt");
    clock_t start = clock();
    pair<Matrix<double>, vector<pair<size_t, size_t>>> LUres = LUdecompos(matrix);
    Matrix<double> LUmatrix = LUres.first;
    cout << "The time: " << (double)(clock() - start) / CLOCKS_PER_SEC << " seconds\n";
    LUmatrix.recordM((char*)"output.txt");

    Matrix<double> b((char*)"inputB.txt");
    Matrix<double> answer = LUsolve(matrix, b);
    answer.recordM((char*)"answer.txt");

    // inverseM(matrix).printM();
    (inverseM(matrix).T() * matrix).printM();


    // Matrix<double> L(matrix.dim), U(matrix.dim);
    // for (size_t i = 0; i < LUmatrix.dim.first; i++) {
    //     for (size_t j = 0; j < LUmatrix.dim.second; j++) {
    //         if (i == j) L[i][j] = 1;
    //         if (i <= j) U[i][j] = LUmatrix[i][j];
    //         else L[i][j] = LUmatrix[i][j];
    //     }
    // }
    // answer = answer.T();
    // cout << '\n'; L.printM(); cout << '\n'; U.printM(); cout << '\n'; (L * U).printM();
    // cout << '\n'; (matrix * answer).T().printM();
    return 0;
}