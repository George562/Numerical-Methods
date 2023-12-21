#include "Jacobi.h"

double epsilon = 0.01;

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    Matrix<double> matrix((char*)"input.txt");
    pair<vector<double>, Matrix<double>> answer = Jacobi(matrix, epsilon);
    answer.second.recordM((char*)"answer.txt");
    for (size_t i = 0; i < matrix.dim.first; i++)
        cout << answer.first[i] << '\n';
    return 0;
}