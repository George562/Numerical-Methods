#include <iostream>
#include "../../Graphica.h"
double epsilon = 0.0001;

float Newton(float (*foo)(float), float (*foo1)(float), float (*foo2)(float), float epsilon, float a, float b) {
    float x = (foo(a) * foo2(a) <= 0 || abs(foo(a) * foo2(a)) >= powf(foo1(a), 2)) ? b : a;

    if (foo(x) * foo2(x) <= 0 || abs(foo(x) * foo2(x)) >= powf(foo1(x), 2)) {
        cerr << "convergence conditions of the method are not met\n"; exit(1);
    }

    float was;
    size_t counter = 0;
    do {
        counter++;
        was = x;
        x = x - foo(x) / foo1(x);
    } while (abs(x - was) > epsilon);
    cout << "counter = " << counter << '\n';
    return x;
}

float Iter(float (*foo)(float), float (*foo1)(float), float epsilon, float a, float b) {
    float q = 0;
    for (float x = a; x <= b; x += (b - a) / 1000) {
        if (foo(x) < a || foo(x) > b) {
            cerr << "convergence conditions of the method are not met\n"; exit(1);
        }
        q = max(q, abs(foo1(x)));
    }
    cout << "q = " << q << '\n';
    float was, x = (b - a) / 2.f;
    size_t counter = 0;
    do {
        counter++;
        was = x;
        x = foo(x);
    } while (abs(x - was) > epsilon);
    cout << "counter = " << counter << '\n';
    return x;
}
