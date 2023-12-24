#include "LinearSolver.h"
#include "../lab1/LU/LU.h"
#include "../lab1/TMA/TMA.hpp"
#include "../lab1/Iter/Iter.h"
#include "../lab1/Jacobi/Jacobi.h"
#include "../lab1/QR/QR.h"
#include <chrono>

int main() {
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    LinearSolver<double> solver;
    Matrix<double> answer;

    char choice; cout << "chouse test: 1 - LU, 2 - TMA, 3 - Iter, 4 - Jacobi, 5 - QR\n"; cout.flush();
    do { cin >> choice; } while (!('1' <= choice || choice <= '5'));

    switch (choice) {
        case '1': { // LU test
            Matrix<double> m1((char*)"LUm.txt");
            Matrix<double> b1((char*)"LUb.txt");

            const auto start1 = std::chrono::high_resolution_clock::now();
            answer = solver.LUsolve(m1, b1);
            const auto end1 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff1 = end1 - start1;
            cout << "The     parallel time: " << diff1.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer1.txt");

            const auto start2 = std::chrono::high_resolution_clock::now();
            answer = LUsolve(m1, b1);
            const auto end2 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff2 = end2 - start2;
            cout << "The non parallel time: " << diff2.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer2.txt");

            break;
        }
        case '2': { // TMA test
            Matrix<double> m2((char*)"TMAm.txt");
            Matrix<double> b2((char*)"TMAb.txt");

            const auto start1 = std::chrono::high_resolution_clock::now();
            answer = solver.TMAsolve(m2, b2);
            const auto end1 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff1 = end1 - start1;
            cout << "The     parallel time: " << diff1.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer1.txt");

            const auto start2 = std::chrono::high_resolution_clock::now();
            answer = TMAsolve(m2, b2);
            const auto end2 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff2 = end2 - start2;
            cout << "The non parallel time: " << diff2.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer2.txt");

            break;
        }
        case '3': { // Iter test
            Matrix<double> m3((char*)"Iterm.txt");
            Matrix<double> b3((char*)"Iterb.txt");
            char TheMethod;
            cout << "What TheMethod do you want to use: 1 - SimpleIter, 2 - Zaydel\n"; cout.flush();
            do { cout << ">>> "; cout.flush(); cin >> TheMethod; } while (TheMethod != '1' && '2' != TheMethod);
            double epsilon; cout << "with epsilon = "; cout.flush(); cin >> epsilon;

            const auto start1 = std::chrono::high_resolution_clock::now();
            answer = (TheMethod == '1') ? solver.SimpleIter(m3, b3, epsilon) : solver.Zaydel(m3, b3, epsilon);
            const auto end1 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff1 = end1 - start1;
            cout << "The     parallel time: " << diff1.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer1.txt");

            const auto start2 = std::chrono::high_resolution_clock::now();
            answer = (TheMethod == '1') ? SimpleIter(m3, b3, epsilon) : Zaydel(m3, b3, epsilon);
            const auto end2 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff2 = end2 - start2;
            cout << "The non parallel time: " << diff2.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer2.txt");

            break;
        }
        case '4': { // Jacobi test
            double epsilon; cout << "with epsilon = "; cout.flush(); cin >> epsilon;
            Matrix<double> m4((char*)"Jacobim.txt");
            
            const auto start1 = std::chrono::high_resolution_clock::now();
            answer = solver.Jacobi(m4, epsilon).second;
            const auto end1 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff1 = end1 - start1;
            cout << "The     parallel time: " << diff1.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer1.txt");

            const auto start2 = std::chrono::high_resolution_clock::now();
            answer = Jacobi(m4, epsilon).second;
            const auto end2 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff2 = end2 - start2;
            cout << "The non parallel time: " << diff2.count() << " seconds\n"; cout.flush();
            answer.recordM((char*)"answer2.txt");

            break;
        }
        case '5': { // QR test
            double epsilon; cout << "with epsilon = "; cout.flush(); cin >> epsilon;
            Matrix<double> m5((char*)"QRm.txt");

            const auto start1 = std::chrono::high_resolution_clock::now();
            auto [R1, C1] = solver.QRsolve(m5, epsilon);
            const auto end1 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff1 = end1 - start1;
            cout << "The     parallel time: " << diff1.count() << " seconds\n"; cout.flush();

            std::ofstream answer1("answer1.txt");
            for (size_t i = 0; i < R1.size(); i++) answer1 << R1[i] << '\n';
            for (size_t i = 0; i < C1.size(); i++)
                answer1 << C1[i].first << " + " << C1[i].second << "i\n";

            const auto start2 = std::chrono::high_resolution_clock::now();
            auto [R2, C2] = QRsolve(m5, epsilon);
            const auto end2 = std::chrono::high_resolution_clock::now();
            const std::chrono::duration<double> diff2 = end2 - start2;
            cout << "The non parallel time: " << diff2.count() << " seconds\n"; cout.flush();

            std::ofstream answer2("answer2.txt");
            for (size_t i = 0; i < R2.size(); i++) answer2 << R2[i] << '\n';
            for (size_t i = 0; i < C2.size(); i++)
                answer2 << C2[i].first << " + " << C2[i].second << "i\n";

            break;
        }
    }

    exit(0);
    return 0;
}