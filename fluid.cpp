#include "fluid.h"
#include "ot.h"

Fluid::Fluid(double volume_fluid, std::vector<Vector> positions) {
    this->volume_fluid = volume_fluid;
    this->positions = positions;
    this->velocities = std::vector<Vector>(positions.size(), Vector(0, 0, 0));
}
void Fluid::simulate_time_step(double dt) {
    int N = positions.size();
    double *lambda = new double[N];
    for (int i = 0; i < N; i++) {
        lambda[i] = volume_fluid / N;
    }
    this->power_diagrams = semidiscrete_ot(positions, lambda, volume_fluid);
    std::cout << "Power diagrams computed" << std::endl;
    std::cout << "Power diagrams size: " << power_diagrams.size() << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << i << " " << power_diagrams[i].area() << std::endl;
        Vector spring =
            (power_diagrams[i].centroid() - positions[i]) / (eps * eps);
        Vector F = spring + mi * Vector(0, -10, 0);
        velocities[i] = velocities[i] + dt * F / mi;
        positions[i] = positions[i] + dt * velocities[i];
    }
    std::cout << "Time step finished" << std::endl;
    delete[] lambda;
}
