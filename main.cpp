#include "polygon.h"
#include "random.h"
#include "svg.h"
#include <cmath>

int main() {
    std::vector<Polygon> show;

    // Polygon t1;
    // t1.vertices.push_back(Vector(0, 0, 0));
    // t1.vertices.push_back(Vector(1, 0, 0));
    // t1.vertices.push_back(Vector(0.5, 1, 0));
    // Polygon t2;
    // t2.vertices.push_back(Vector(1, 1, 0));
    // t2.vertices.push_back(Vector(0, 1, 0));
    // t2.vertices.push_back(Vector(0.5, 0, 0));
    // show.push_back(Polygon::sutherland_hodgman(t1, t2));
    Polygon rand_pol;
    for (int i = 0; i < 1000; i++) {
        rand_pol.vertices.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }

    Polygon square;
    square.vertices.push_back(Vector(0, 0, 0));
    square.vertices.push_back(Vector(1, 0, 0));
    square.vertices.push_back(Vector(1, 1, 0));
    square.vertices.push_back(Vector(0, 1, 0));

    show = Polygon::venoroi(rand_pol, square);
    save_svg(show, "image.svg");
    return 0;
}
