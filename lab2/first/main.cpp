#include "tools.h"

float foo(float x) { return sin(x) - x * x + 1.f; }
float foo1(float x) { return cos(x) - 2 * x; }
float foo2(float x) { return - sin(x) - 2.f; }

int main() {
    makeGraph([](double x){ return sin(x) + 1; });
    makeGraph([](double x){ return x * x; });
    sf::Thread thread(ShowGraphics); thread.launch();
    float a, b; cout << "a = "; cin >> a; cout << "b = "; cin >> b;
    cout << "Choose method: 1 - Newton, 2 - Iter\n> ";
    int choose; cin >> choose;
    float res;
    if (choose == 1)
        cout << Newton(foo, foo1, foo2, epsilon, a, b) << '\n';
    else
        cout << Iter([](float x) { return powf(sin(x) + 1, 0.5); },
                     [](float x) { return 0.5f * cos(x) / powf(sin(x) + 1, 0.5); },
                     epsilon, a, b) << '\n';
    thread.wait();
    return 0;
}