#include "polygon.h"
class Fluid {
    void simulate_time_step(double dt);
    std::vector<Vector> positions;
    std::vector<Vector> velocities;
    double volume_fluid;
    double eps = 0.004;
    double dt = 0.002;
    double mi = 200;

  public:
    Fluid(double volume_fluid, std::vector<Vector> velocities);
};
