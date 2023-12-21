#include "tools.h"
#include "../../Graphica.h"

int main() {
    ifstream input("input.txt");
    size_t n; input >> n;
    vector<double> X(n), foo(n);
    for (size_t i = 0; i < n; i++) input >> X[i];
    for (size_t i = 0; i < n; i++) input >> foo[i];
    for (size_t i = 0; i < n; i++) addPoint(X[i], foo[i]);
    SetDiapazon(X[0], X[n - 1]);
    MaxY = 6; MinY = 0;

    sf::Thread thread(ShowGraphics); thread.launch();
    sf::Mutex mutex;
    int dim;
    while (true) {
        cout << "input dim: "; cin >> dim;
        if (dim == 0) break;
        double (*res)(double) = Approximation(foo, X, dim);
        double Sum_of_Squared_Errors = 0;
        for (size_t i = 0; i < n; i++) Sum_of_Squared_Errors += pow(foo[i] - res(X[i]), 2);
        cout << "Sum of Squared Errors = " << Sum_of_Squared_Errors << '\n';
        mutex.lock();
        makeGraph(res, sf::Color(rand() % 256, rand() % 256, rand() % 256));
        mutex.unlock();
    }
    thread.wait();
    return 0;
}

/*
Интерполяция — нахождение неизвестных промежуточных значений некоторой функции,
по имеющемуся дискретному набору её известных значений

Сплайн — это кусочно заданная функция, то есть совокупность нескольких функций,
каждая из которых задана на каком-то множестве значений аргумента, причём эти множества попарно непересекающиеся.

Аппроксимация - замена одних объектов другими, в каком-то смысле близкими к исходным, но более простыми.
*/