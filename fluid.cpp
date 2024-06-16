#include "fluid.h"
#include "ot.h"

Fluid::Fluid(double volume_fluid, std::vector<Vector> positions) {
    this->volume_fluid = volume_fluid;
    this->positions = positions;
    this->velocities = std::vector<Vector>(positions.size(), Vector(0, 0, 0));
}
void Fluid::simulate_time_step(double dt, double *weights) {
    int N = positions.size();
    double *lambda = new double[N + 1];
    for (int i = 0; i < N; i++) {
        lambda[i] = volume_fluid / N;
    }
    lambda[N] = 1 - volume_fluid;
    this->power_diagrams = semidiscrete_ot(positions, lambda, volume_fluid, weights);
    for (int i = 0; i < N; i++) {
        Vector spring =
            (power_diagrams[i].centroid() - positions[i]) / (eps * eps);
        Vector F = spring + mi * Vector(0, -10, 0);
        velocities[i] = velocities[i] + dt * F / mi;
        positions[i] = positions[i] + dt * velocities[i];
        if (positions[i][0] < 0) {
            positions[i][0] = 0;
            velocities[i][0] = 0;
        }
        if (positions[i][1] < 0) {
            positions[i][1] = 0;
            velocities[i][1] = 0;
        }
        if (positions[i][0] > 1) {
            positions[i][0] = 1;
            velocities[i][0] = 0;
        }
        if (positions[i][1] > 1) {
            positions[i][1] = 1;
            velocities[i][1] = 0;
        }
    }
    delete[] lambda;
}
