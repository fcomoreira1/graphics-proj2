#include "random.h"

double uniform_distribution() {
    static thread_local std::mt19937 generator;
    std::uniform_real_distribution<double> distribution(0, 1);
    return distribution(generator);
}

Vector random_direction() {
    double r1 = 2 * M_PI * uniform_distribution(), r2 = uniform_distribution();
    double auxr2 = sqrt(r2 * (1 - r2));

    return Vector(cos(r1) * auxr2,
                  sin(r1) * auxr2,
                  1 - 2 * r2);
}

// Code from Wikipedia on BoxMuller
std::pair<double, double> box_muller(double mu, double sigma) {
    constexpr double two_pi = 2.0 * M_PI;

    // create two random numbers, make sure u1 is greater than zero
    double u1 = uniform_distribution(), u2 = uniform_distribution();

    // compute z0 and z1
    auto mag = sigma * sqrt(-2.0 * log(u1));
    auto z0 = mag * cos(two_pi * u2) + mu;
    auto z1 = mag * sin(two_pi * u2) + mu;

    return std::make_pair(z0, z1);
}
