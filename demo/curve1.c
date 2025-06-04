/*
 * Library:  lmfit (Levenberg-Marquardt least squares fitting)
 *
 * File:     demo/curve1.c
 *
 * Contents: Example for curve fitting with lmcurve():
 *           fit a data set y(x) by a curve f(x;p).
 *
 * Note:     Any modification of this example should be copied to
 *           the manual page source lmcurve.pod and to the wiki.
 *
 * Author:   Joachim Wuttke <j.wuttke@fz-juelich.de> 2004-2013
 *
 * Licence:  see ../COPYING (FreeBSD)
 *
 * Homepage: apps.jcns.fz-juelich.de/lmfit
 */

#include "lmcurve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* model function: a parabola */

double f(const double t, const double *p)
{
    return p[0] + p[1] * t + p[2] * t * t;
}

void run_test(const char *test_name, lm_control_struct *control, double *par,
              int n, int m, double *t, double *y)
{
    printf("\n========================================\n");
    printf("TEST: %s\n", test_name);
    printf("========================================\n");

    /* Reset parameters to bad starting values */
    par[0] = 100.0;
    par[1] = 0.0;
    par[2] = -10.0;

    lm_status_struct status;

    /* Call lmcurve */
    lmcurve(n, par, m, t, y, f, control, &status);

    printf("\nResults:\n");
    printf("status after %d function evaluations:\n  %s\n",
           status.nfev, lm_infmsg[status.outcome]);

    printf("obtained parameters:\n");
    for (int i = 0; i < n; ++i)
        printf("  par[%i] = %12.5f\n", i, par[i]);
    printf("obtained norm:\n  %12.6g\n", status.fnorm);

    printf("fitting data as follows:\n");
    for (int i = 0; i < m; ++i)
        printf("  t[%2d]=%4g y=%6g fit=%10.5f residue=%12.6g\n",
               i, t[i], y[i], f(t[i], par), y[i] - f(t[i], par));

    /* Expected solution */
    double expected[3] = {0.12987, -0.05, 1.03052};
    double tolerance = 1e-3;      /* Tolerance for parameter accuracy */
    double norm_tolerance = 1e-3; /* Tolerance for norm accuracy */

    /* Check parameter accuracy */
    int params_correct = 1;
    printf("\nParameter accuracy check:\n");
    for (int i = 0; i < n; i++)
    {
        double error = fabs(par[i] - expected[i]);
        printf("  par[%d]: obtained = %.5f, expected = %.5f, error = %.2e",
               i, par[i], expected[i], error);
        if (error > tolerance)
        {
            printf(" [FAIL]");
            params_correct = 0;
        }
        else
        {
            printf(" [OK]");
        }
        printf("\n");
    }

    /* Check norm accuracy */
    double expected_norm = 0.450685;
    double norm_error = fabs(status.fnorm - expected_norm);
    printf("Norm check: obtained=%.6f, expected=%.6f, error=%.2e",
           status.fnorm, expected_norm, norm_error);
    int norm_correct = (norm_error < norm_tolerance);
    if (norm_correct)
    {
        printf(" [OK]\n");
    }
    else
    {
        printf(" [FAIL]\n");
    }

    /* Overall success criteria */
    int convergence_success = (status.outcome <= 4);         /* Convergence achieved */
    int accuracy_success = (params_correct && norm_correct); /* Correct solution */

    printf("\nOverall assessment:\n");
    printf("  Convergence: %s (outcome %d)\n",
           convergence_success ? "SUCCESS" : "FAILURE", status.outcome);
    printf("  Accuracy: %s\n", accuracy_success ? "SUCCESS" : "FAILURE");

    if (convergence_success && accuracy_success)
    {
        printf("OVERALL: SUCCESS\n");
    }
    else if (convergence_success && !accuracy_success)
    {
        printf("OVERALL: CONVERGED TO WRONG SOLUTION\n");
    }
    else if (!convergence_success && accuracy_success)
    {
        printf("OVERALL: CORRECT SOLUTION BUT POOR CONVERGENCE\n");
    }
    else
    {
        printf("OVERALL: FAILURE\n");
    }
}

int main()
{
    int n = 3;     /* number of parameters in model function f */
    double par[3]; /* parameters array */

    /* data points: a slightly distorted standard parabola */
    int m = 9;
    double t[9] = {-4., -3., -2., -1., 0., 1., 2., 3., 4.};
    double y[9] = {16.6, 9.9, 4.4, 1.1, 0., 1.1, 4.2, 9.3, 16.4};

    /* Expected solution from polynomial regression:
     * Y = A + B1*X + B2*X^2
     * A  = 0.12987 ± 0.09299
     * B1 = -0.05   ± 0.02375
     * B2 = 1.03052 ± 0.01048
     *
     * Correct answer: p = {0.12987, -0.05, 1.03052}
     */

    printf("Expected solution: [0.12987, -0.05, 1.03052]\n");

    /* Base control structure */
    lm_control_struct control = lm_control_double;
    control.verbosity = 3;
    control.ftol = 1e-8;
    control.xtol = 1e-8;
    control.gtol = 1e-8;
    control.patience = 1000; /* Allow enough iterations */

    /* TEST 1: No scaling (scale_diag = 0) */
    control.scale_diag = 0;
    control.bounds = NULL;
    run_test("No Scaling (scale_diag = 0)", &control, par, n, m, t, y);

    /* TEST 2: Auto scaling (scale_diag = 1) */
    control.scale_diag = 1;
    control.bounds = NULL;
    run_test("Auto Scaling (scale_diag = 1)", &control, par, n, m, t, y);

    /* TEST 3: User-provided scales only */
    control.scale_diag = 0; /* Will be ignored when user scales provided */

    /* Create bounds structure with user-provided scales */
    static lm_bounds_struct bounds3;
    bounds3.lower = NULL;
    bounds3.upper = NULL;
    bounds3.bound_type = NULL;

    /* User-provided scales: different scaling for each parameter */
    static double user_scales[3] = {0.1, 0.01, 1.0};
    bounds3.scales = user_scales;

    control.bounds = &bounds3;
    run_test("User-Provided Scales [0.1, 0.01, 1.0]", &control, par, n, m, t, y);

    /* TEST 4: Bounds only (no user scales, with auto scaling) */
    control.scale_diag = 1; /* Enable auto scaling with bounds */

    static lm_bounds_struct bounds4;
    static double lower_bounds4[3] = {-1.0, -0.5, 0.5};
    static double upper_bounds4[3] = {5.0, 0.5, 2.0};
    static int bound_types4[3] = {LM_BOUND_BOTH, LM_BOUND_BOTH, LM_BOUND_BOTH};

    bounds4.lower = lower_bounds4;
    bounds4.upper = upper_bounds4;
    bounds4.bound_type = bound_types4;
    bounds4.scales = NULL; /* No user scales, allow auto scaling */

    control.bounds = &bounds4;
    run_test("Bounds Only + Auto Scaling", &control, par, n, m, t, y);

    /* TEST 5: Bounds + User-provided scales */
    control.scale_diag = 1; /* Will be ignored when user scales provided */

    static lm_bounds_struct bounds5;
    static double lower_bounds5[3] = {-1.0, -0.5, 0.5};
    static double upper_bounds5[3] = {5.0, 0.5, 2.0};
    static int bound_types5[3] = {LM_BOUND_BOTH, LM_BOUND_BOTH, LM_BOUND_BOTH};
    static double user_scales5[3] = {2.0, 0.5, 1.5}; /* Different scales */

    bounds5.lower = lower_bounds5;
    bounds5.upper = upper_bounds5;
    bounds5.bound_type = bound_types5;
    bounds5.scales = user_scales5;

    control.bounds = &bounds5;
    run_test("Bounds + User-Provided Scales [2.0, 0.5, 1.5]", &control, par, n, m, t, y);

    /* TEST 6: Lower bounds only (all parameters have lower bounds) */
    control.scale_diag = 0;

    static lm_bounds_struct bounds6;
    static double lower_bounds6[3] = {-1000.0, -1000.0, -1000.0};
    static double upper_bounds6[3] = {5.0, 0.5, 2.0}; /* Not used for lower bounds only */
    static int bound_types6[3] = {LM_BOUND_LOWER, LM_BOUND_LOWER, LM_BOUND_LOWER};

    bounds6.lower = lower_bounds6;
    bounds6.upper = upper_bounds6;
    bounds6.bound_type = bound_types6; /* Will be ignored */
    bounds6.scales = NULL;             /* No scaling */

    control.bounds = &bounds6;
    run_test("Lower Bounds Only: all parameters have lower bounds [-1000, +inf)", &control, par, n, m, t, y);

    printf("\n========================================\n");
    printf("ALL TESTS COMPLETED\n");
    printf("========================================\n");

    return 0;
}