#pragma once
#include "lbfgs/lbfgs.h"
#include "vector.h"
#include <iostream>
#include <vector>

const bool ENABLE_KDTREE = false;
const int NCirc = 1000;

class Polygon {
    std::vector<Vector> vertices;

  public:
    Polygon() : vertices() {}
    Polygon(std::vector<Vector> vertices) : vertices(vertices) {}
    Polygon clip(Vector u, Vector v) const;
    double area() const;
    Vector centroid() const;
    Vector operator[](int i) const {
        if (i < 0 or i >= vertices.size()) {
            std::cerr << "Polygon index out of range" << std::endl;
            exit(1);
        }
        return vertices[i];
    }
    Vector &operator[](int i) {
        if (i < 0 or i >= vertices.size()) {
            std::cerr << "Polygon index out of range" << std::endl;
            exit(1);
        }
        return vertices[i];
    }
    int size() const { return vertices.size(); }
    static double area_triangle(Vector A, Vector B, Vector C);
};
Polygon sutherland_hodgman(const Polygon &subject, const Polygon &clipper);
std::vector<Polygon> venoroi(const std::vector<Vector> &points,
                             const Polygon &default_shape = Polygon(std::vector(
                                 {Vector(0, 0, 0), Vector(1, 0, 0),
                                  Vector(1, 1, 0), Vector(0, 1, 0)})));
std::vector<Polygon> power_diagrams(
    const std::vector<Vector> &points, const lbfgsfloatval_t *weights,
    const bool use_air = false,
    const Polygon &default_shape = Polygon(std::vector(
        {Vector(0, 0, 0), Vector(1, 0, 0), Vector(1, 1, 0), Vector(0, 1, 0)})));
