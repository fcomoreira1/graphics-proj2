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
    save_svg(show, "images/image_sutherland_hodgman.svg");
}

void test_vonoroi(int num_points = 2000) {
    std::cout << "Testing Vonoroi" << std::endl;
    std::vector<Vector> random_points;
    for (int i = 0; i < num_points; i++) {
        random_points.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    std::vector<Polygon> show = venoroi(random_points);
    save_svg(show, random_points, "images/image_vonoroi.svg");
    std::cout << "Image saved in images/image_vonoroi.svg" << std::endl;
}

void test_optim_power_diagrams(int num_points = 200){
    std::vector<Vector> random_points;
    for (int i = 0; i < num_points; i++) {
        random_points.push_back(
            Vector(uniform_distribution(), uniform_distribution(), 0));
    }
    double *lambda = new double[num_points];
    Vector C(0.5, 0.5, 0);
    double s;
    for (int i = 0; i < num_points; i++) {
        lambda[i] = exp(-(random_points[i] - C).norm2() / 0.02);
        // lambda[i] = 1.0 / num_points;
        s += lambda[i];
    }
    for (int i = 0; i < random_points.size(); i++) {
        lambda[i] /= s;
    }
    std::vector<Polygon> show = semidiscrete_ot(random_points, lambda);
    delete[] lambda;
    save_svg(show, random_points, "images/image_power_diagrams.svg");
}

void test_fluid() {
    int NPoints = 150;
    int num_frames = 600;
    double expected_volume = 0.4;
    double dt = 0.003;

    std::vector<Vector> random_points(NPoints);
    for (int i = 0; i < NPoints; i++) {
        random_points[i] =
            Vector(uniform_distribution(), uniform_distribution(), 0);
    }
    std::cerr << "Creating fluid" << std::endl;
    Fluid fluid(expected_volume, random_points);
    double *weights = new double[NPoints + 1];
    for (int i = 0; i < NPoints; i++) {
        weights[i] = 1.0;
    }
    weights[NPoints] = 0.0;
    for (int i = 0; i < num_frames; i++) {
        fluid.simulate_time_step(dt, weights);
        save_svg_animated(fluid.power_diagrams, "images/image_fluid.svg", i, num_frames);
        // save_frame(fluid.power_diagrams, "image_fluid", i);
        std::cerr << "Saved frame " << i << std::endl;
    }
    free(weights);
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        std::cerr << "Usage: ./main [0: sutherland_hodgman | 1 number_of_points: vonoroi | 2 number_of_points: "
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
        if (argc == 2){
            std::cerr << "number of points missing" << std::endl;
            exit(-1);
        }
        test_vonoroi(std::stoi(argv[2]));
        break;
    case 2:
        if (argc == 2){
            std::cerr << "number of points missing" << std::endl;
            exit(-1);
        }
        test_optim_power_diagrams(std::stoi(argv[2]));
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
