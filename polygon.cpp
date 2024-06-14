#include "polygon.h"
#include "kdtree.h"
#include "lbfgs/lbfgs.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <omp.h>

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
    if (this->size() == 0)
        return Vector(0, 0, 0);
    if (this->size() == 1)
        return vertices[0];
    if (this->size() == 2)
        return (vertices[0] + vertices[1]) / 2;
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
        Vector u = clipper[i];
        Vector v = clipper[(i + 1) % clipper.size()];
        Vector N = Vector((v - u).data[1], -(v - u).data[0], 0);
        N.normalize();
        output = output.clip(u, N);
    }
    return output;
}

std::vector<Polygon> power_diagrams(std::vector<Vector> &points,
                                    const lbfgsfloatval_t *weights,
                                    const bool use_air,
                                    const Polygon &default_shape) {
    int N = points.size();
    std::vector<Polygon> Pow(N);
    KDTree kdtree(points);
    auto clip_cell = [&points, &weights](int i, int j, Polygon &cell) {
        if (i == j)
            return;
        Vector M = 0.5 * (points[i] + points[j]);
        Vector Norm = points[j] - points[i];
        Vector Mprime =
            M + (weights[i] - weights[j]) * (Norm) / (2 * Norm.norm2());
        Norm.normalize();
        cell = cell.clip(Mprime, Norm);
    };
#pragma omp parallel
    for (int i = 0; i < N; i++) {
        Polygon cell = default_shape;

        if (ENABLE_KDTREE) {
            std::vector<u_long> ret_indexes;
            int k = init_k;
            while (true) {
                kdtree.findNeighbors(points[i].data, k, ret_indexes);
                for (auto j : ret_indexes)
                    clip_cell(i, j, cell);
                double max_dist_cell = -1;
                for (int j = 0; j < cell.size(); j++) {
                    max_dist_cell =
                        std::max(max_dist_cell, (points[i] - cell[j]).norm());
                }
                double cur_dist =
                    (points[ret_indexes[ret_indexes.size() - 1]] - points[i])
                        .norm();
                if (k >= N || cur_dist > 2 * max_dist_cell)
                    break;
                k = 2 * k;
            }
        } else {
            for (int j = 0; j < N; j++)
                clip_cell(i, j, cell);
        }

        if (use_air) {
            if (weights[N] > weights[i]) {
                cell = Polygon();
            } else {
                double radius = std::sqrt(weights[i] - weights[N]);
                for (int j = 0; j < NCirc; j++) {
                    double ang1 = 2 * M_PI * j / NCirc;
                    double ang2 = 2 * M_PI * (j + 1) / NCirc;
                    Vector u =
                        points[i] + radius * Vector(cos(ang1), sin(ang1), 0);
                    Vector v =
                        points[i] + radius * Vector(cos(ang2), sin(ang2), 0);
                    Vector Norm = (u + v) / 2 - points[i];
                    Norm.normalize();
                    cell = cell.clip(u, Norm);
                }
            }
        }
#pragma omp critical
        { Pow[i] = cell; }
    }
    return Pow;
}

std::vector<Polygon> venoroi(std::vector<Vector> &points,
                             const Polygon &default_shape) {
    lbfgsfloatval_t *weights = new lbfgsfloatval_t[points.size()];
    memset(weights, 0, points.size() * sizeof(lbfgsfloatval_t));
    auto ret = power_diagrams(points, weights, false, default_shape);
    delete[] weights;
    return ret;
}
