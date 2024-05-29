#pragma once
#include "polygon.h"

struct OptimizationInstance {
    std::vector<Vector> points;
    double *lambda;
};

std::vector<Polygon> optimal_transport(const std::vector<Vector> &points,
                                       double *lambda);
