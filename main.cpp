#include "polygon.h"
#include "random.h"
#include "svg.h"
#include <cmath>

void test_sutherland_hodgman() {
    std::vector<Polygon> show;

    Polygon t1;
    t1.vertices.push_back(Vector(0, 0, 0));
    t1.vertices.push_back(Vector(1, 0, 0));
    t1.vertices.push_back(Vector(0.5, 1, 0));
    Polygon t2;
    t2.vertices.push_back(Vector(1, 1, 0));
    t2.vertices.push_back(Vector(0, 1, 0));
    t2.vertices.push_back(Vector(0.5, 0, 0));
    show.push_back(sutherland_hodgman(t1, t2));
    save_svg(show, "image_sutherland_hodgman.svg");
}

void test_vonoroi() {
    Polygon rand_pol;
    for (int i = 0; i < 2000; i++) {
        rand_pol.vertices.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    std::vector<Polygon> show = venoroi(rand_pol);
    save_svg(show, "image_vonoroi.svg");
}

void test_power_diagrams() {
    Polygon rand_pol;
    for (int i = 0; i < 1000; i++) {
        rand_pol.vertices.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    double *lambda = new double[rand_pol.vertices.size()];
    Vector C(0.5, 0.5, 0);
    for (int i = 0; i < rand_pol.vertices.size(); i++) {
        lambda[i] = exp(-(rand_pol.vertices[i] - C).norm2() / 0.02);
    }
    std::vector<Polygon> show = optim_power_diagrams(rand_pol, lambda);
    delete[] lambda;
    save_svg(show, "image_power_diagrams.svg");
}

int main() {
    test_power_diagrams();
    return 0;
}
