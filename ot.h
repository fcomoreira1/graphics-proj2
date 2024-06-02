#pragma once
#include "polygon.h"

struct OptimizationInstance {
    const std::vector<Vector> points;
    const double *lambda;
    const double desired_volume_fluid;
};

std::vector<Polygon> semidiscrete_ot(const std::vector<Vector> &points,
                                     const double *lambda,
                                     const double desired_volume = 1.0);
