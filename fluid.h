#include "polygon.h"
class Fluid {
    std::vector<Vector> velocities;
    double volume_fluid;
    double eps = 0.004;
    double mi = 200;

  public:
    std::vector<Vector> positions;
    std::vector<Polygon> power_diagrams;
    void simulate_time_step(double dt, double *weights);
    Fluid(double volume_fluid, std::vector<Vector> positions);
};
