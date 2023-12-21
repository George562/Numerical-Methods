#include <omp.h>
#include <chrono>
#include <iostream>
// #include <thread>

int main() {
    size_t res = 1;
    const auto start1 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for reduction(+:res)
    for (size_t i = 1; i < 1000000000; i++) {
        res = i + res;
    }
    const auto end1 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> diff1 = end1 - start1;
    std::cout << "The time: " << diff1.count() << " seconds\n";
    std::cout << res << '\n';
    res = 1;
    const auto start2 = std::chrono::high_resolution_clock::now();
    for (size_t i = 1; i < 1000000000; i++) {
        res = i + res;
    }
    const auto end2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "The time: " << diff2.count() << " seconds\n";
    std::cout << res << '\n';
    return 0;
}