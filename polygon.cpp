#include "polygon.h"
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

double Polygon::area() const {
    double sum = 0;
    for (int i = 0; i < vertices.size(); i++) {
        Vector A = vertices[i];
        Vector B = vertices[(i + 1) % vertices.size()];
        sum += A[0] * B[1] - A[1] * B[0];
    }
    return 0.5 * sum;
}

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

vector<Polygon> power_diagrams(const Polygon &P, const lbfgsfloatval_t *weights,
                               const Polygon &default_shape) {
    std::vector<Polygon> Pow(P.size());
#pragma omp parallel
    for (int i = 0; i < P.size(); i++) {
        Polygon cell = default_shape;
        for (int j = 0; j < P.size(); j++) {
            if (i == j)
                continue;
            Vector M = 0.5 * (P[i] + P[j]);
            Vector N = P[j] - P[i];
            Vector Mprime =
                M + (weights[i] - weights[j]) * (N) / (2 * dot(N, N));
            N.normalize();
            cell = cell.clip(M, N);
        }
#pragma omp critical
        { Pow[i] = cell; }
    }
    return Pow;
}

std::vector<Polygon> venoroi(const Polygon &P, const Polygon &default_shape) {
    double *weights = new lbfgsfloatval_t[P.size()];
    memset(weights, 0, P.size() * sizeof(lbfgsfloatval_t));
    auto ret = power_diagrams(P, weights, default_shape);
    delete[] weights;
    return ret;
}

struct PolygonOptimizationInstance {
    Polygon P;
    double *lambda;
};

/**
 * @brief Computes the integral of ||P_i - x||^2 over a polygon.
 */
static double integrate_cell(const Polygon &P, const Vector &X) {
    double par_integral = 0;
    std::vector<Vector> c(3);
    c[0] = P[0];
    for (int j = 1; j < P.size() - 1; j++) {
        c[1] = P[j];
        c[2] = P[j + 1];
        double par_sum = 0.0;
        for (int k = 0; k < 3; k++) {
            for (int l = k; l < 3; l++) {
                par_sum += dot(c[k] - X, c[l] - X);
            }
        }
        double area_triangle = 0.5 * cross(c[1] - c[0], c[2] - c[0]).norm();
        par_integral += std::abs(par_sum * area_triangle / 6.0);
    }
    return par_integral;
}

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *weights,
                                lbfgsfloatval_t *grad_g, const int n,
                                const lbfgsfloatval_t step) {

    PolygonOptimizationInstance *inst = (PolygonOptimizationInstance *)instance;
    std::vector<Polygon> Pow = power_diagrams(inst->P, weights);
    double g = 0.0;

    for (int i = 0; i < inst->P.size(); i++) {
        // Flipped signs because lbfgs minimizes the function
        g -= integrate_cell(Pow[i], inst->P.vertices[i]) -
             std::abs(Pow[i].area()) * weights[i] +
             inst->lambda[i] * weights[i];
        grad_g[i] = (lbfgsfloatval_t)std::abs(Pow[i].area()) - inst->lambda[i];
    }
    return (lbfgsfloatval_t)g;
}

static int progress(void *instance, const lbfgsfloatval_t *x,
                    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
                    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
                    const lbfgsfloatval_t step, int n, int k, int ls) {
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

static int semi_disc_optimal_transport(const Polygon &P, double *lambda,
                                       lbfgsfloatval_t *weights) {
    lbfgs_parameter_t param;
    lbfgsfloatval_t fx;
    PolygonOptimizationInstance instance = {P, lambda};
    lbfgs_parameter_init(&param);
    int ret =
        lbfgs(P.size(), weights, &fx, evaluate, progress, &instance, &param);
    if (ret != 0) {
        std::cerr << "LBFGS optimization failed with code: " << ret
                  << std::endl;
    }
    return ret;
}

std::vector<Polygon> optim_power_diagrams(const Polygon &P, double *lambda) {
    lbfgsfloatval_t *weights = lbfgs_malloc(P.size());
    if (weights == nullptr) {
        std::cerr << "ERROR: Failed to allocate a memory block for variables."
                  << std::endl;
        exit(1);
    }
    semi_disc_optimal_transport(P, lambda, weights);
    auto ret = power_diagrams(P, weights);
    lbfgs_free(weights);
    return ret;
}
