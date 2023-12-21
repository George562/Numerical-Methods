// Генератор матриц размера n x n
// при запуске генератора с консоли введите через пробел размер матрицы
// Пример: ./generator.exe 2
// генерация чисел происходит в диапазоне от 1 до max
// max моно ввести третьей переменной. если не вводить, то max = RAND_MAX
#include <iostream>
#include <string>
#include <ctime>
#include <vector>
using namespace std;

int main(int arg, char *kwargs[]) {
    srand(time(0)); rand();
    size_t n = stoull(kwargs[1]), max = stoull(kwargs[2]);
    max = (max == 0) ? RAND_MAX : max;
    cout << n << ' ' << n << '\n';
    std::vector<std::vector<size_t>> v(n, std::vector<size_t>(n));

    for (size_t i = 0; i < n; i++)
        for (size_t j = i; j < n; j++)
            v[i][j] = v[j][i] = 1 + rand() % max;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++)
            cout << v[i][j] << '\t';
        cout << '\n';
    }
    return 0;
}