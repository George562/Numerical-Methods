// Генератор матриц размера n x m
// при запуске генератора с консоли введите через пробел размер матрицы
// Пример: ./generator.exe 2 2
// генерация чисел происходит в диапазоне от 1 до max
// max моно ввести третьей переменной. если не вводить, то max = RAND_MAX
#include <iostream>
#include <string>
#include <ctime>
using namespace std;

int main(int arg, char *kwargs[]) {
    srand(time(0)); rand();
    size_t n = stoull(kwargs[1]), m = stoull(kwargs[2]), max = stoull(kwargs[3]);
    max = (max == 0) ? RAND_MAX : max;
    cout << n << ' ' << m << '\n';
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++)
            cout << 1 + rand() % max << '\t';
        cout << '\n';
    }
    return 0;
}