#include "ot.h"
#include "polygon.h"
#include "random.h"
#include "svg.h"
#include <chrono>
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
    std::vector<Vector> random_poitns;
    for (int i = 0; i < 2000; i++) {
        random_poitns.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    std::vector<Polygon> show = venoroi(random_poitns);
    save_svg(show, "image_vonoroi.svg");
}

void test_power_diagrams() {
    int test_size = 1000;

    std::vector<Vector> vertices;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < test_size; i++) {
        Vector temp =
            Vector(distribution(generator), distribution(generator), 0);
        vertices.push_back(temp);
    }

    // Polygon edges;
    // edges.vertices.push_back(Vector(0., 1., 0.));
    // edges.vertices.push_back(Vector(1., 1., 0.));
    // edges.vertices.push_back(Vector(1., 0., 0.));
    // edges.vertices.push_back(Vector(0., 0., 0.));

    double *weights = new double[test_size];
    for (int i = 0; i < test_size; i++) {
        weights[i] = (distribution(generator) / 4.0);
    }

    // vertices.push_back(Vector(0., 1., 0.));
    // vertices.push_back(Vector(1., 1., 0.));
    // vertices.push_back(Vector(1., 0., 0.));
    // vertices.push_back(Vector(0., 0., 0.));

    auto start_time = std::chrono::system_clock::now();
    save_svg(venoroi(vertices), vertices, "test1_unweight.svg");
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> run_time = end_time - start_time;
    std::cout << "Finished rendering 1!" << std::endl
              << "Rendering time: " << run_time.count() << " s" << std::endl
              << std::endl;

    start_time = std::chrono::system_clock::now();
    save_svg(power_diagrams(vertices, weights), vertices, "test1_weight.svg");
    end_time = std::chrono::system_clock::now();
    run_time = end_time - start_time;
    std::cout << "Finished rendering 1 weighted!" << std::endl
              << "Rendering time: " << run_time.count() << " s" << std::endl
              << std::endl;
}

void test_optim_power_diagrams() {
    std::vector<Vector> random_points;
    for (int i = 0; i < 1000; i++) {
        random_points.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    double *lambda = new double[random_points.size()];
    Vector C(0.5, 0.5, 0);
    double s;
    for (int i = 0; i < random_points.size(); i++) {
        lambda[i] = exp(-(random_points[i] - C).norm2() / 0.02);
        s += lambda[i];
    }
    for (int i = 0; i < random_points.size(); i++) {
        lambda[i] /= s;
    }
    std::vector<Polygon> show = optimal_transport(random_points, lambda);
    delete[] lambda;
    save_svg(show, random_points, "image_power_diagrams.svg");
}

int main() {
    auto start_time = std::chrono::system_clock::now();
    test_optim_power_diagrams();
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> run_time = end_time - start_time;
    std::cout << "Rendering time: " << run_time.count() << " s" << std::endl
              << std::endl;
    return 0;
}
