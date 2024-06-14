#pragma once
#include <random>
#include "vector.h"

double uniform_distribution();
std::pair<double, double> box_muller(double mu, double sigma);
Vector random_direction();
