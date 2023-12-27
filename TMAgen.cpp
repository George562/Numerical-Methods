// Генератор матриц размера 3 x n
// при запуске генератора с консоли введите размер матрицы
// Пример: ./generator.exe 2
// генерация чисел происходит в диапазоне от 1 до max
// max можно ввести второй переменной. если не вводить, то max = RAND_MAX
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
using namespace std;

int main(int arg, char *kwargs[]) {
    srand(time(0)); rand();
    size_t n = stoull(kwargs[1]), max = stoull(kwargs[2]);
    max = (max == 0) ? RAND_MAX : max;
    cout << "3 " << n << '\n';
    std::vector<std::vector<size_t>> v(3, std::vector<size_t>(n));

    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < n; j++) {
            v[i][j] = 1 + rand() % max;
        }
    
    cout << "0\t";
    for (size_t i = 1; i < n; i++) cout << v[0][i] << '\t';
    cout << '\n';
    for (size_t i = 0; i < n; i++) cout << v[1][i] + v[0][i] + v[2][i] + rand() % (max / 100) << '\t';
    cout << '\n';
    for (size_t i = 0; i < n - 1; i++) cout << v[2][i] << '\t';
    cout << "0\n";
    return 0;
}