#include "fluid.h"
#include "ot.h"
#include "polygon.h"
#include "random.h"
#include "svg.h"
#include <chrono>
#include <cmath>

void test_sutherland_hodgman() {
    std::vector<Polygon> show;
    Polygon t1(
        std::vector{Vector(0, 0, 0), Vector(1, 0, 0), Vector(0.5, 1, 0)});
    Polygon t2(
        std::vector{Vector(1, 1, 0), Vector(0, 1, 0), Vector(0.5, 0, 0)});
    show.push_back(sutherland_hodgman(t1, t2));
    save_svg(show, "image_sutherland_hodgman.svg");
}

void test_vonoroi() {
    std::cout << "Testing Vonoroi" << std::endl;
    std::vector<Vector> random_points;
    for (int i = 0; i < 2000; i++) {
        random_points.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    std::vector<Polygon> show = venoroi(random_points);
    save_svg(show, "image_vonoroi.svg");
    std::cout << "Image saved in image_vonoroi.svg" << std::endl;
}

void test_optim_power_diagrams() {
    std::vector<Vector> random_points;
    for (int i = 0; i < 200; i++) {
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
    std::vector<Polygon> show = semidiscrete_ot(random_points, lambda);
    delete[] lambda;
    save_svg(show, random_points, "image_power_diagrams.svg");
}

void test_fluid() {
    int NPoints = 20;
    double expected_volume = 0.02;
    double dt = 0.02;
    std::vector<Vector> random_points(NPoints);
    for (int i = 0; i < NPoints; i++) {
        random_points[i] =
            Vector(uniform_distribution(), uniform_distribution(), 0);
    }
    std::cerr << "Creating fluid" << std::endl;
    Fluid fluid(expected_volume, random_points);
    for (int i = 0; i < 100; i++) {
        fluid.simulate_time_step(dt);
        save_svg(fluid.power_diagrams, fluid.positions,
                 "fluid/image_fluid" + std::to_string(i) + ".svg");
        // save_frame(fluid.power_diagrams, "image_fluid", i);
        std::cerr << "Saved frame " << i << std::endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        std::cerr << "Usage: ./main [0: sutherland_hodgman | 1: vonoroi | 2: "
                     "optim_power_diagrams | 3: fluid]"
                  << std::endl;
        exit(1);
    }
    auto start_time = std::chrono::system_clock::now();
    switch (std::stoi(argv[1])) {
    case 0:
        test_sutherland_hodgman();
        break;
    case 1:
        test_vonoroi();
        break;
    case 2:
        test_optim_power_diagrams();
        break;
    case 3:
        test_fluid();
        break;
    }
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> run_time = end_time - start_time;
    std::cout << "Rendering time: " << run_time.count() << " s" << std::endl
              << std::endl;
    return 0;
}
