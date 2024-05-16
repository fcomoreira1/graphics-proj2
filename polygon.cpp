#include "polygon.h"
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

Polygon Polygon::sutherland_hodgman(const Polygon &subject,
                                    const Polygon &clipper) {
    Polygon output = subject;
    for (int i = 0; i < clipper.vertices.size(); i++) {
        Polygon input = output;
        output.vertices.clear();
        Vector u = clipper.vertices[i];
        Vector v = clipper.vertices[(i + 1) % clipper.vertices.size()];
        Vector N = Vector((v - u).data[1], -(v - u).data[0], 0);
        N.normalize();
        output = output.clip(u, N);
    }
    return output;
}

std::vector<Polygon> Polygon::venoroi(const Polygon &P,
                                      const Polygon &default_shape) {
    std::vector<Polygon> veno(P.vertices.size());
#pragma omp parallel
    for (int i = 0; i < P.vertices.size(); i++) {
        Polygon cell = default_shape;
        for (int j = 0; j < P.vertices.size(); j++) {
            if (i == j)
                continue;
            Vector M = 0.5 * (P.vertices[i] + P.vertices[j]);
            Vector N = P.vertices[j] - P.vertices[i];
            N.normalize();
            cell = cell.clip(M, N);
        }
#pragma omp critical
        { veno[i] = cell; }
    }
    return veno;
}
