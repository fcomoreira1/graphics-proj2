#pragma once
#include "polygon.h"

struct OptimizationInstance {
    std::vector<Vector> points;
    const double *lambda;
    const double desired_volume_fluid;
};

std::vector<Polygon> semidiscrete_ot(std::vector<Vector> &points,
                                     const double *lambda,
                                     const double desired_volume = 1.0,
                                     double *weights = nullptr);
