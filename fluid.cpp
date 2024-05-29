#include "fluid.h"
#include "ot.h"

Fluid::Fluid(double volume_fluid, std::vector<Vector> velocities) {
    this->volume_fluid = volume_fluid;
    this->velocities = velocities;
    this->positions = std::vector<Vector>(velocities.size(), Vector(0, 0, 0));
}
void Fluid::simulate_time_step(double dt) {
    double *lambda = new double[positions.size()];
    for (int i = 0; i < positions.size(); i++) {
        lambda[i] = volume_fluid / positions.size();
    }
    std::vector<Polygon> Pow = optimal_transport(positions, lambda);
    for (int i = 0; i < positions.size(); i++) {
        Vector spring = (Pow[i].centroid() - positions[i]) / (eps * eps);
        Vector F = spring + mi * Vector(0, -10, 0);
        velocities[i] = velocities[i] + dt * F / mi;
        positions[i] = positions[i] + dt * velocities[i];
    }
}
