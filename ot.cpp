#include "ot.h"
#include "random.h"
/**
 * @brief Computes the integral of ||P_i - x||^2 over a polygon.
 */
static double integrate_cell(const Polygon &P, const Vector &X) {
    if (P.size() < 3)
        return 0.0;
    double par_integral = 0;
    for (int j = 1; j < P.size() - 1; j++) {
        Vector A = P[0] - X, B = P[j] - X, C = P[j + 1] - X;
        double par_sum = dot(A, B) + dot(B, C) + dot(A, C) + dot(A, A) +
                         dot(B, B) + dot(C, C);

        par_integral += std::abs(
            par_sum * Polygon::area_triangle(P[0], P[j], P[j + 1]) / 6.0);
    }
    return par_integral;
}

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *weights,
                                lbfgsfloatval_t *grad_g, const int n,
                                const lbfgsfloatval_t step) {

    OptimizationInstance *inst = (OptimizationInstance *)instance;
    std::vector<Polygon> Pow = power_diagrams(inst->points, weights);
    double g = 0.0;
    for (int i = 0; i < n; i++) {
        // Flipped signs because lbfgs minimizes the function
        grad_g[i] =
            (lbfgsfloatval_t)(std::abs(Pow[i].area()) - inst->lambda[i]);
        g -= integrate_cell(Pow[i], inst->points[i]) - grad_g[i] * weights[i];
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

std::vector<Polygon> optimal_transport(const std::vector<Vector> &points,
                                       double *lambda) {
    lbfgsfloatval_t *weights = lbfgs_malloc(points.size());
    for (int i = 0; i < points.size(); i++) {
        weights[i] = uniform_distribution() / 10;
    }
    lbfgs_parameter_t param;
    lbfgsfloatval_t fx;
    OptimizationInstance instance = {points, lambda};
    lbfgs_parameter_init(&param);
    param.max_iterations = 1000;
    int ret = lbfgs(points.size(), weights, &fx, evaluate, progress,
                    (void *)&instance, &param);
    if (ret < 0) {
        std::cerr << "LBFGS optimization failed with code: " << ret
                  << std::endl;
    }
    // semi_disc_optimal_transport(points, lambda, weights);
    auto Pow = power_diagrams(points, weights);
    lbfgs_free(weights);
    return Pow;
}
