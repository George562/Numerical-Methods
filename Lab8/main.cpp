#include "tools.h"
#include <thread>
#include "../Graphica.h"
#include "../UI/bar.h"

double f         (double x, double y, double t) { return 0.                                    ; }

double leftFunc  (double x, double y, double t) { return              sinh(y) * exp(-3 * a * t); }
double rightFunc (double x, double y, double t) { return -2 *         sinh(y) * exp(-3 * a * t); }

double downFunc  (double x, double y, double t) { return cos(2 * x) *           exp(-3 * a * t); }
double upFunc    (double x, double y, double t) { return cos(2 * x) * 0.75 *    exp(-3 * a * t); }

double initFunc  (double x, double y, double t) { return cos(2 * x) * sinh(y)                  ; }

double answer    (double x, double y, double t) { return cos(2 * x) * sinh(y) * exp(-3 * a * t); }

int main() {
    // параболический тип
    ios::sync_with_stdio(false); cin.tie(0); cout.tie(0);
    setlocale(LC_ALL, "rus");

    setA(1.);

    loadFonts();

    Bar<int> barY;
    barY.ValueText.setFont(font);
    barY.setColors(sf::Color(0, 0, 0), sf::Color(250, 120, 0), sf::Color(120, 120, 120));
    barY.setWidth(scw);
    barY.setPosition(0, sch - STANDART_BAR_HEIGHT * 2);
    drawableStuff.push_back(&barY);

    Bar<int> barT;
    barT.ValueText.setFont(font);
    barT.setColors(sf::Color(0, 0, 0), sf::Color(250, 120, 0), sf::Color(120, 120, 120));
    barT.setWidth(scw);
    barT.setPosition(0, sch - STANDART_BAR_HEIGHT);
    drawableStuff.push_back(&barT);

    PlacedText textY, textT;
    textY.setString("y = "); textT.setString("t = ");
    textY.setPosition(barY.ValueText.getPosition() - sf::Vector2f(textY.getSize().x, 0.f));
    textT.setPosition(barT.ValueText.getPosition() - sf::Vector2f(textT.getSize().x, 0.f));
    drawableStuff.push_back(&textY); drawableStuff.push_back(&textT);

    thread thread1(ShowGraphics);

    double Lx = M_PI_4, Ly = log(2.), Lt = 0.5;
    BC leftCondition    {0., 1., leftFunc };
    BC rightCondition   {1., 0., rightFunc};
    BC downCondition    {1., 0., downFunc };
    BC upCondition      {0., 1., upFunc   };
    BC initCondition    {0., 1., initFunc };

    cout << "какой метод использовать?\n1 - переменных направлений\n2 - дробных шагов\n>"; cout.flush();
    int ans; cin >> ans;

    auto resFunc = ans == 1 ? AlternatingDirectionMethod : FractionalStepsMethod;

    cout << "с максимальным разбиением Nx Ny Nt = "; cout.flush();
    size_t Nx, Ny, Nt; cin >> Nx >> Ny >> Nt;
    vector<vector<vector<double>>> res = resFunc(Nx, Ny, Nt, Lx, Ly, Lt, leftCondition, rightCondition, downCondition, upCondition, initCondition, f);

    vector<vector<vector<double>>> answerM(Nt + 1, vector<vector<double>>(Nx + 1, vector<double>(Ny + 1)));
    double error = 0.;
    for (size_t i = 0; i < Nt + 1; i++)
        for (size_t j = 0; j < Nx + 1; j++)
            for (size_t k = 0; k < Ny + 1; k++) {
                answerM[i][j][k] = answer((j * Lx) / Nx, (k * Ly) / Ny, (i * Lt) / Nt);
                error = max(error, abs(answerM[i][j][k] - res[i][j][k]));
            }

    cout << "error = " << error << '\n'; cout.flush();

    vector<pair<double, double>> g1(Nx + 1);
    for (int i = 0; i < g1.size(); i++)
        g1[i] = {(i * Lx) / Nx, res[Nt][i][Ny]};

    vector<pair<double, double>> g2(5001);
    for (int i = 0; i < g2.size(); i++) {
        g2[i].first = (i * Lx) / 5000;
        g2[i].second = answer((i * Lx) / 5000, Ly, Lt);
    }
    makeGraphByPoints(g1);
    makeGraphByPoints(g2, sf::Color::Green);

    Scale<int> scaleT{0, (int)Nt, (int)Nt};
    barT.setValue(scaleT);
    Scale<int> scaleY{0, (int)Ny, (int)Ny};
    barY.setValue(scaleY);

    thread thread2([&]{
        cout << "---------------------------------------------------------------\n"
             << "|  N  :     hx      :     hy      :     tau     :    error    |\n"
             << "+-----+-------------+-------------+-------------+-------------+\n";

        vector<pair<double, double>> g3(50 - 3 + 1);
        double error;
        vector<vector<vector<double>>> res;
        for (int n = 3; n <= 50; n++) {
            error = 0;
            res = resFunc(n, n, n, Lx, Ly, Lt, leftCondition, rightCondition, downCondition, upCondition, initCondition, f);

            for (size_t i = 0; i < n + 1; i++)
                for (size_t j = 0; j < n + 1; j++)
                    for (size_t k = 0; k < n + 1; k++)
                        error = max(error, abs(answer((j * Lx) / n, (k * Ly) / n, (i * Lt) / n) - res[i][j][k]));
            g3[n - 3] = {Lx / n, error};

            cout << "| " << setw(3) << n << " : " << setw(11) << Lx / n << " : " << setw(11)
                 << Ly / n << " : " << setw(11) << Lt / n << " : " << setw(11) << error << " |\n"
                 << "+-----+-------------+-------------+-------------+-------------+\n";
        }
        cout << '\n'; cout.flush();
        makeGraphByPoints(g3);
    });
    thread thread3([&]{
        while (window == nullptr || window->isOpen()) {
            if (window->hasFocus() && sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                if (sf::Mouse::getPosition(*window).y >= barT.getPosition().y || sf::Mouse::getPosition(*window).y >= barY.getPosition().y) {
                    MouseBuffer = Mouse.getPosition(*window);
                    if (sf::Mouse::getPosition(*window).y >= barT.getPosition().y) {
                        scaleT.cur = (double)scaleT.top * (double)sf::Mouse::getPosition(*window).x / (double)scw;
                        normalize(scaleT);
                    } else {
                        scaleY.cur = (double)scaleY.top * (double)sf::Mouse::getPosition(*window).x / (double)scw;
                        normalize(scaleY);
                    }
                    textY.setPosition(barY.ValueText.getPosition() - sf::Vector2f(textY.getSize().x, 0.f));
                    textT.setPosition(barT.ValueText.getPosition() - sf::Vector2f(textT.getSize().x, 0.f));

                    for (int i = 0; i < g1.size(); i++)
                        graph[0][i].y = -res[scaleT.cur][i][scaleY.cur];

                    for (int i = 0; i < g2.size(); i++)
                        graph[1][i].y = -answer((i * Lx) / 5000, Ly * (double)scaleY.cur / (double)scaleY.top, Lt * (double)scaleT.cur / (double)scaleT.top);
                }
            }
        }
    });
    thread1.join();
    thread2.join();
    thread3.join();
    return 0;
}