#include "tools.h"
#include <thread>
#include "../UI/bar.h"

double f         (double x, double y) { return 0.        ; }

double leftFunc  (double x, double y) { return 0.        ; }
double rightFunc (double x, double y) { return y         ; }

double downFunc  (double x, double y) { return sin(x)    ; }
double upFunc    (double x, double y) { return 0.        ; }

double answer    (double x, double y) { return y * sin(x); }

int main() {
    // элиптический тип
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    setlocale(LC_ALL, "rus");

    setC(1.);

    loadFonts();
    Bar<int> bar;
    bar.ValueText.setFont(font);
    bar.setColors(sf::Color(0, 0, 0), sf::Color(250, 120, 0), sf::Color(120, 120, 120));
    bar.setWidth(scw);
    bar.setPosition(0, sch - STANDART_BAR_HEIGHT);
    drawableStuff.push_back(&bar);
    thread thread1(ShowGraphics);

    double Lx = M_PI_2, Ly = 1.;
    BC leftCondition  {0.,  1., leftFunc };
    BC rightCondition {0.,  1., rightFunc};
    BC downCondition  {1.,  0., downFunc };
    BC upCondition    {1., -1., upFunc   };

    cout << "какой метод использовать?\n1 - Либмана\n2 - Зейделя\n3 - Верхней релаксации\n>"; cout.flush();
    int ans; cin >> ans;
    IterType::IterType type = IterType::IterType(ans - 1);
    cout << "epsilon = "; cout.flush();
    double epsilon; cin >> epsilon;

    double omega = 0.;
    if (ans == 3) { cout << "omega = "; cout.flush(); cin >> omega; }

    cout << "с максимальным разбиением Nx Ny = "; cout.flush();
    size_t Nx, Ny; cin >> Nx >> Ny;
    Matrix<double> res = Iter(Nx, Ny, Lx, Ly, epsilon, leftCondition, rightCondition, downCondition, upCondition, f, type, omega);
    size_t K = res.dim.second - 1;

    vector<pair<double, double>> g1(Nx + 1);
    for (int i = 0; i < g1.size(); i++)
        g1[i] = {(i * Lx) / Nx, res[i][K]};

    vector<pair<double, double>> g2(5001);
    for (int i = 0; i < g2.size(); i++) {
        g2[i].first = (i * Lx) / 5000;
        g2[i].second = answer((i * Lx) / 5000, Ly);
    }
    makeGraphByPoints(g1);
    makeGraphByPoints(g2, sf::Color::Green);

    Scale<int> scale{0, (int)K, (int)K};
    bar.setValue(scale);

    // thread thread2([&]{
    //     cout << "----------------------------------------------------------\n"
    //          << "|  Nx  :   Ny  :      hx     :     hy      :    error    |\n"
    //          << "+-----+--------+-------------+-------------+-------------+\n";

    //     vector<pair<double, double>> g3(Nx - 3 + 1);
    //     double error;
    //     Matrix<double> res;
    //     for (int n = 3; n <= min(Nx, Ny); n++) {
    //         error = 0;
    //         res = Iter(n, n, Lx, Ly, epsilon / n, leftCondition, rightCondition, downCondition, upCondition, f, type, omega);

    //         for (int i = 0; i <= n; i++)
    //             for (int j = 0; j <= n; j++)
    //                 error += max(abs(answer((i * Lx) / n, (j * Ly) / n) - res[i][j]) - error, 0.);
    //         g3[n - 3] = {Lx / n, error};

    //         cout << "| " << setw(3) << n << " : " << setw(6) << n
    //              << " : " << setw(11) << Lx / n << " : " << setw(11) << Ly / n
    //              << " : " << setw(11) << error << " |\n"
    //              << "+-----+--------+-------------+-------------+-------------+\n";
    //     }
    //     cout << '\n'; cout.flush();
    //     makeGraphByPoints(g3);
    // });
    thread thread3([&]{
        while (window == nullptr || window->isOpen()) {
            if (window->hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                if (sf::Mouse::getPosition(*window).y >= bar.getPosition().y) {
                    MouseBuffer = Mouse.getPosition(*window);
                    scale.cur = (double)scale.top * (double)sf::Mouse::getPosition(*window).x / (double)scw;
                    normalize(scale);

                    for (int i = 0; i < g1.size(); i++)
                        graph[0][i].y = -res[i][scale.cur];

                    for (int i = 0; i < g2.size(); i++)
                        graph[1][i].y = -answer((i * Lx) / 5000, Ly * (double)scale.cur / (double)scale.top);
                }
            }
        }
    });
    thread1.join();
    // thread2.join();
    thread3.join();
    return 0;
}