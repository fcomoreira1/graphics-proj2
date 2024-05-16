#ifndef POLYGON_H
#define POLYGON_H
#include "vector.h"
#include <vector>

class Polygon {
  public:
    std::vector<Vector> vertices;
    Polygon() : vertices() {}
    static Polygon sutherland_hodgman(const Polygon &subject,
                                      const Polygon &clipper);
    Polygon clip(Vector u, Vector v) const;
    static std::vector<Polygon> venoroi(const Polygon &P,
                                        const Polygon &default_shape);
};
#endif
