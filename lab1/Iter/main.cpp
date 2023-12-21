#include "Iter.h"

double epsilon = 0.01;

int main() {
    Matrix<double> matrix((char*)"inputM.txt");
    Matrix<double> b((char*)"inputB.txt");
    char TheMethod;
    cout << "What TheMethod do you want to use?\n"
    "1 - SimpleIter\n"
    "2 - Zaydel\n";
    do { cout << ">>> "; cin >> TheMethod; } while (TheMethod != '1' && '2' != TheMethod);
    Matrix<double> answer;
    if (TheMethod == '1') answer = SimpleIter(matrix, b, epsilon);
    else                  answer = Zaydel(matrix, b, epsilon);
    answer.recordM((char*)"answer.txt");
    (matrix * answer.T()).printM(); // проверка
    return 0;
}