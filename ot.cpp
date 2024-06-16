#include "ot.h"
#include "lbfgs/lbfgs.h"
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
                                lbfgsfloatval_t *grad, const int n,
                                const lbfgsfloatval_t step) {

    OptimizationInstance *inst = (OptimizationInstance *)instance;
    std::vector<Polygon> Pow =
        power_diagrams(inst->points, weights, inst->desired_volume_fluid < 1.0);
    double g = 0.0;
    double est_vol_air = 1.0;
    for (int i = 0; i < n - 1; i++) {
        // Flipped signs because lbfgs minimizes the function
        est_vol_air -= std::abs(Pow[i].area());
        grad[i] = std::abs(Pow[i].area()) - inst->lambda[i];
        g -= integrate_cell(Pow[i], inst->points[i]) - grad[i] * weights[i];
    }
    double des_volume_air = 1 - inst->desired_volume_fluid;
    g -= weights[n - 1] * (des_volume_air - est_vol_air);
    grad[n - 1] = -(des_volume_air - est_vol_air);
    return g;
}

static int progress(void *instance, const lbfgsfloatval_t *x,
                    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
                    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
                    const lbfgsfloatval_t step, int n, int k, int ls) {
    static double prev_fx;
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    if (k > 200 && (prev_fx - fx) < 1e-6) {
        prev_fx = -1;
        return 1;
    }
    prev_fx = fx;
    return 0;
}

std::vector<Polygon> semidiscrete_ot(std::vector<Vector> &points,
                                     const double *lambda,
                                     const double desired_volume) {
    const int N = points.size();
    lbfgsfloatval_t *weights = lbfgs_malloc(N + 1);
    for (int i = 0; i < N; i++)
        // weights[i] = uniform_distribution() / 10;
        weights[i] = 1.0;
    weights[N] = 0.0;

    lbfgsfloatval_t fx;
    OptimizationInstance instance = {points, lambda, desired_volume};

    lbfgs_parameter_t param;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
    param.max_iterations = 1000;
    lbfgs_parameter_init(&param);

    int ret = lbfgs(N + 1, weights, &fx, evaluate, progress, (void *)&instance,
                    &param);
    if (ret < 0 && ret != LBFGSERR_MAXIMUMITERATION) {
        std::cerr << "LBFGS optimization failed with code: " << ret
                  << std::endl;
    }
    auto Pow = power_diagrams(points, weights, desired_volume < 1.0);
    std::cout << "Dumping weights" << std::endl;
    double total_weight = 0.0;
    for (int i = 0; i < N; i++) {
        std::cout << "area_" << i << " -> "
                  << Pow[i].area() << std::endl;
        total_weight += Pow[i].area();
    }
    std::cout << "Area air -> " << 1 - total_weight << std::endl;

    lbfgs_free(weights);
    return Pow;
}
