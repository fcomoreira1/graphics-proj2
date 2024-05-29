#include "polygon.h"
#include "lbfgs/lbfgs.h"
#include "random.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <omp.h>

using std::vector;

Polygon Polygon::clip(Vector u, Vector N) const {
    Polygon output;
    for (int j = 0; j < vertices.size(); j++) {
        Vector B = vertices[j];
        Vector A = vertices[j > 0 ? j - 1 : vertices.size() - 1];
        double t = dot(u - A, N) / dot(B - A, N);
        bool Binside = dot(B - u, N) <= 0;
        bool Ainside = dot(A - u, N) <= 0;

        if (Binside) { // B inside
            if (!Ainside) {
                output.vertices.push_back(A + (B - A) * t);
            }
            output.vertices.push_back(B);
        } else if (Ainside) {
            output.vertices.push_back(A + (B - A) * t);
        }
    }
    return output;
}

double Polygon::area_triangle(Vector A, Vector B, Vector C) {
    return 0.5 * cross(B - A, C - A).norm();
}

double Polygon::area() const {
    if (vertices.size() < 3)
        return 0.0;
    double sum = 0;
    for (int i = 1; i < vertices.size() - 1; i++) {
        sum += area_triangle(vertices[0], vertices[i], vertices[i + 1]);
    }
    return sum;
}
Vector Polygon::centroid() const {
    double a = area();
    Vector c;
    for (int i = 0; i < this->size(); i++) {
        int ip = (i + 1) % this->size();
        double xy =
            vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1];
        c = c + Vector((vertices[i][0] + vertices[ip][0]) * xy,
                       (vertices[i][1] + vertices[ip][1]) * xy, 0);
    }
    return c / (6 * a);
};

Polygon sutherland_hodgman(const Polygon &subject, const Polygon &clipper) {
    Polygon output = subject;
    for (int i = 0; i < clipper.size(); i++) {
        Polygon input = output;
        output.vertices.clear();
        Vector u = clipper[i];
        Vector v = clipper[(i + 1) % clipper.size()];
        Vector N = Vector((v - u).data[1], -(v - u).data[0], 0);
        N.normalize();
        output = output.clip(u, N);
    }
    return output;
}

vector<Polygon> power_diagrams(const std::vector<Vector> &points,
                               const lbfgsfloatval_t *weights,
                               const Polygon &default_shape) {
    std::vector<Polygon> Pow(points.size());
#pragma omp parallel
    for (int i = 0; i < points.size(); i++) {
        Polygon cell = default_shape;
        for (int j = 0; j < points.size(); j++) {
            if (i == j)
                continue;
            Vector M = 0.5 * (points[i] + points[j]);
            Vector N = points[j] - points[i];
            Vector Mprime =
                M + (weights[i] - weights[j]) * (N) / (2 * N.norm2());
            N.normalize();
            cell = cell.clip(Mprime, N);
        }
        // #pragma omp critical
        // { Pow[i] = cell; }
        std::swap(Pow[i], cell);
        // Pow[i] = std::move(cell);
    }
    return Pow;
}

std::vector<Polygon> venoroi(const std::vector<Vector> &points,
                             const Polygon &default_shape) {
    lbfgsfloatval_t *weights = new lbfgsfloatval_t[points.size()];
    memset(weights, 0, points.size() * sizeof(lbfgsfloatval_t));
    auto ret = power_diagrams(points, weights, default_shape);
    delete[] weights;
    return ret;
}
