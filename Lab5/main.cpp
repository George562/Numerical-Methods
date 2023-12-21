#include "tools.h"
#include <thread>
#include "../UI/bar.h"

double firstFunc  (double x, double t) { return          sin(t); }
double secondFunc (double x, double t) { return         -sin(t); }
double initialFunc(double x, double t) { return              0.; }
double answer     (double x, double t) { return sin(t) * cos(x); }
double f          (double x, double t) { return cos(x) * (sin(t) + cos(t)); }

int main() {
    // параболический тип
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    setlocale(LC_ALL, "Russian");

    thread thread1(ShowGraphics);

    loadFonts();
    Bar<int> bar;
    bar.ValueText.setFont(font);
    bar.setColors(sf::Color(0, 0, 0), sf::Color(250, 120, 0), sf::Color(120, 120, 120));
    bar.setWidth(scw);
    bar.setPosition(0, sch - STANDART_BAR_HEIGHT);
    drawableStuff.push_back(&bar);

    double L = M_PI_2, T = 0.6;
    BC leftCondition     {0., 1., firstFunc  };
    BC rightCondition    {1., 0., secondFunc };
    BC initialConditions {0., 1., initialFunc};

    cout << "какой метод использовать?\n1 - явный\n2 - неявный\n3 - смешанный\n>"; cout.flush();
    int ans; cin >> ans;
    double theta = ans == 1 ? 0. : ans == 2 ? 1. : 0.;
    cout << "с максимальным разибением N = "; cout.flush();
    size_t N; cin >> N;
    Matrix<double> res = finiteDifferenceMethod(N, L, T, leftCondition, rightCondition, initialConditions, f, theta);
    size_t K = res.dim.second - 1;

    vector<pair<double, double>> g1(N + 1);
    for (int i = 0; i < g1.size(); i++)
        g1[i] = {(i * L) / N, res[i][K]};

    vector<pair<double, double>> g2(5001);
    for (int i = 0; i < g2.size(); i++) {
        g2[i].first = (i * L) / 5000;
        g2[i].second = answer((i * L) / 5000, T);
    }
    makeGraphByPoints(g1);
    makeGraphByPoints(g2, sf::Color::Green);

    Scale<int> scale{0, (int)K, (int)K};
    bar.setValue(scale);

    thread thread2([&]{
        cout << "------------------------------------------------------------------\n"
             << "|  N  :   K    :      h      :     tau     : sigma :    error    |\n"
             << "+-----+--------+-------------+-------------+-------+-------------+\n";

        vector<pair<double, double>> g3(N - 3 + 1);
        double error;
        Matrix<double> res;
        for (int n = 3; n <= N; n++) {
            error = 0;
            res = finiteDifferenceMethod(n, L, T, leftCondition, rightCondition, initialConditions, f, theta);
            size_t K = res.dim.second - 1;

            for (int i = 0; i <= n; i++)
                for (int j = 0; j <= K; j++)
                    error += max(abs(answer((i * L) / n, (j * T) / K) - res[i][j]) - error, 0.);
            g3[n - 3] = {L / n, error};

            cout << "| " << setw(3) << n << " : " << setw(6) << size_t(T / (0.495 * pow(L / n, 2.)))
                 << " : " << setw(11) << L / n << " : " << setw(11) << T / (T / (0.495 * pow(L / n, 2.)))
                 << " : " << setw(5) << (T / (T / (0.495 * pow(L / n, 2.)))) / pow(L / n, 2.) << " : "
                 << setw(11) << error << " |\n"
                 << "+-----+--------+-------------+-------------+-------+-------------+\n";
        }
        cout << '\n'; cout.flush();
        makeGraphByPoints(g3);
    });
    thread thread3([&]{
        while (window == nullptr || window->isOpen()) {
            if (window->hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                if (sf::Mouse::getPosition(*window).y >= bar.getPosition().y) {
                    MouseBuffer = Mouse.getPosition(*window);
                    scale.cur = (double)scale.top * (double)sf::Mouse::getPosition(*window).x / (double)scw;
                    normalize(scale);

                    res = finiteDifferenceMethod(N, L, T, leftCondition, rightCondition, initialConditions, f, theta);
                    for (int i = 0; i < g1.size(); i++)
                        graph[0][i].y = -res[i][scale.cur];

                    vector<pair<double, double>> g2(5001);
                    for (int i = 0; i < g2.size(); i++) {
                        graph[1][i].y = -answer((i * L) / 5000, T * (double)scale.cur / (double)scale.top);
                    }
                }
            }
        }
    });
    thread1.join();
    thread2.join();
    thread3.join();
    return 0;
}