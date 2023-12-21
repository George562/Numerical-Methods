#include <cmath>
#include "SFML-2.5.1/include/SFML/Graphics.hpp"
// -lsfml-graphics -lsfml-window -lsfml-system
using namespace std;
using vp = vector<sf::Vector2f>;
using vvp = vector<vp>;
sf::RenderWindow* window = nullptr;
int scw = 800, sch = 800;
sf::Event event;
sf::Vector2f CameraPos(scw / 2, sch / 2);
sf::Vector2i MouseBuffer;
sf::Mouse Mouse;
sf::Font font; sf::Text AxisText;
sf::CircleShape PointShape;
float pointRad = 5.f;
int countOfStep = 5000;
float Start = -3.25, End = 3.25;
float step = (End - Start) / countOfStep;
float MaxX = End + 0.25, MinX = Start - 0.25;
vvp graph(0);
vp points(0);
vector<sf::Color> colors(0);
float zoom = pow(1.1, 50);
bool showPoints = true;
stringstream stringStream;
float stepBy = 1;
float stepsBy[3] = {0.5, 0.4, 0.5};
int IndexStepBy = 0;
sf::VertexArray VertexOfAxisLine(sf::Lines, 2), line(sf::Lines, 2);
vector<sf::Drawable*> drawableStuff;


float ConvertX(float x) { return (x - CameraPos.x) / zoom; }
float ConvertY(float y) { return (y - CameraPos.y) / zoom; }
float XAxisPlace(sf::Text text) {
    if (ConvertY(0) < 0 && ConvertY(sch - text.getGlobalBounds().height - 15) > 0)
        return 6.f;
    return ((ConvertY(0) >= 0) ? 0 : sch - text.getGlobalBounds().height - 10) - CameraPos.y;
}
// ------------------------------------------------------------
float YAxisPlace(sf::Text text) {
    if (ConvertX(-7) < 0 && ConvertX(scw - text.getGlobalBounds().width - 7) > 0)
        return 6.f;
    return ((ConvertX(-7) >= 0) ? 0 : scw - text.getGlobalBounds().width - 2) - CameraPos.x;
}
// ------------------------------------------------------------
void SetAxisTextByFloatX(float x) {
    stringStream.str(""); stringStream << x;
    AxisText.setString(stringStream.str());
    AxisText.setPosition(sf::Vector2f(x * zoom, XAxisPlace(AxisText)) + CameraPos);
}
// ------------------------------------------------------------
void SetAxisTextByFloatY(float y) {
    stringStream.str(""); stringStream << y;
    AxisText.setString(stringStream.str());
    AxisText.setPosition(sf::Vector2f(YAxisPlace(AxisText), -y * zoom) + CameraPos);
}
// ------------------------------------------------------------
void SetVertexOfAxisLineX(float x) {
    VertexOfAxisLine[0] = sf::Vertex(sf::Vector2f(CameraPos.x + x * zoom, 0), sf::Color::Black);
    VertexOfAxisLine[1] = sf::Vertex(sf::Vector2f(CameraPos.x + x * zoom, sch), sf::Color::Black);
}
// ------------------------------------------------------------
void SetVertexOfAxisLineY(float y) {
    VertexOfAxisLine[0] = sf::Vertex(sf::Vector2f(0, CameraPos.y - y * zoom), sf::Color::Black);
    VertexOfAxisLine[1] = sf::Vertex(sf::Vector2f(scw, CameraPos.y - y * zoom), sf::Color::Black);
}

void makeGraph(double (*foo)(double), sf::Color color=sf::Color::Blue) {
    colors.push_back(color);
    graph.push_back({});
    for (double x = Start; x <= End; x += step)
        graph[graph.size() - 1].push_back({float(x), float(-foo(x))});
}
void makeGraphByPoints(vector<pair<double, double>> foo, sf::Color color=sf::Color::Blue) {
    colors.push_back(color);
    graph.push_back({});
    for (pair<double, double>& x: foo)
        graph[graph.size() - 1].push_back({float(x.first), float(-x.second)});
}

void addPoint(double x, double y) {
    points.emplace_back(x, -y);
}

void ShowGraphics() {
    window = new sf::RenderWindow(sf::VideoMode(scw, sch), "Graphick");
    font.loadFromFile("arial.ttf");
    AxisText.setFont(font);
    AxisText.setCharacterSize(14);
    AxisText.setFillColor(sf::Color::Black);

    PointShape.setFillColor(sf::Color::Green);
    PointShape.setRadius(pointRad);
    while (window->isOpen()) {
        if (window->hasFocus()) {
            if (Mouse.isButtonPressed(Mouse.Left))
                CameraPos += sf::Vector2f(Mouse.getPosition(*window) - MouseBuffer);
        }
        MouseBuffer = Mouse.getPosition(*window);
        window->clear(sf::Color::White);
        AxisText.setPosition(0, 0);
        stringStream.str(""); stringStream <<  "x: " << (MouseBuffer.x - CameraPos.x) / zoom <<
                                             "\ny: " << (-MouseBuffer.y + CameraPos.y) / zoom;
        AxisText.setString(stringStream.str());
        window->draw(AxisText);
        // axis
        for (float x = min(-1, int(ConvertX(scw) / stepBy)) * stepBy; x >= ConvertX(0); x -= stepBy) { // x < 0
            SetAxisTextByFloatX(x); SetVertexOfAxisLineX(x);
            window->draw(AxisText); window->draw(VertexOfAxisLine);
        }
        for (float x = max(0, int(ConvertX(0) / stepBy)) * stepBy; x <= ConvertX(scw); x += stepBy) { // x >= 0
            SetAxisTextByFloatX(x); SetVertexOfAxisLineX(x);
            window->draw(AxisText); window->draw(VertexOfAxisLine);
        }
        for (float y = min(-1, int(-ConvertY(0) / stepBy)) * stepBy; y >= -ConvertY(sch); y -= stepBy) { // y < 0
            SetAxisTextByFloatY(y); SetVertexOfAxisLineY(y);
            window->draw(AxisText); window->draw(VertexOfAxisLine);
        }
        for (float y = max(1, int(-ConvertY(sch) / stepBy)) * stepBy; y <= -ConvertY(0); y += stepBy) { // y > 0
            SetAxisTextByFloatY(y); SetVertexOfAxisLineY(y);
            window->draw(AxisText); window->draw(VertexOfAxisLine);
        }
        SetVertexOfAxisLineY(0);
        window->draw(VertexOfAxisLine);
        // graph
        for (int j = 0; j < graph.size(); j++)
            for (int i = 0; i < graph[j].size() - 1; i++) {
                line[0] = sf::Vertex(sf::Vector2f(graph[j][i].x, graph[j][i].y) * zoom + CameraPos, colors[j]);
                line[1] = sf::Vertex(sf::Vector2f(graph[j][i + 1].x, graph[j][i + 1].y) * zoom + CameraPos, colors[j]);
                window->draw(line);
            }
        // points
        if (showPoints)
            for (int i = 0; i < points.size(); i++) {
                PointShape.setPosition(sf::Vector2f{points[i].x - pointRad / zoom,
                    points[i].y - pointRad / zoom} * zoom + CameraPos);
                window->draw(PointShape);
            }

        while (window->pollEvent(event))
            if (event.type == sf::Event::Closed) window->close();
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Escape) window->close();
                if (event.key.code == sf::Keyboard::H) showPoints = !showPoints;
            } else if (event.type == sf::Event::MouseWheelScrolled) {
                float coef = pow(1.1, event.mouseWheelScroll.delta);
                zoom *= coef;
                CameraPos = (sf::Vector2f)MouseBuffer - ((sf::Vector2f)MouseBuffer - CameraPos) * coef;

                if (stepBy * zoom > 190) {
                    stepBy = stepBy * stepsBy[IndexStepBy];
                    IndexStepBy = (IndexStepBy + 1) % 3;
                } else if (stepBy * zoom < 80) {
                    IndexStepBy = (IndexStepBy + 2) % 3;
                    stepBy /= stepsBy[IndexStepBy];
                }
            }
        for (int i = 0; i < drawableStuff.size(); i++) {
            window->draw(*drawableStuff[i]);
        }
        window->display();
    }
}

void SetDiapazon(float s, float e) {
    Start = s, End = e;
    step = (End - Start) / countOfStep;
    MaxX = ceil(End) + 0.25, MinX = floor(Start) - 0.25;
}