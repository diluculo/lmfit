#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lmmin.h"

#define N 70
#define PI 3.14159265358979323846

/* --- Complex number structure --- */
typedef struct
{
    double real;
    double imag;
} Complex;

/* --- Simple Randles circuit parameters --- */
typedef struct
{
    double Rs, Rct, Cdl;
} Params;

/* --- Global data arrays --- */
static double frequencies[N];
static Complex impedances[N];

/* Weighting method enumeration */
typedef enum
{
    WEIGHT_UNITY = 0,         /* No weighting (unit weights) */
    WEIGHT_MODULUS_MEAS,      /* Weight by 1/|Z_measured| */
    WEIGHT_PROPORTIONAL_MEAS, /* Weight by 1/|Re|, 1/|Im| separately */
    WEIGHT_MODULUS_FIT,       /* Weight by 1/|Z_fitted| */
    WEIGHT_PROPORTIONAL_FIT,  /* Weight by 1/|Re_fit|, 1/|Im_fit| separately */
} WeightingMethod;

/* Global weighting state */
static WeightingMethod current_weighting = WEIGHT_UNITY; /* Default to no weighting */

/* --- Model functions --- */
Complex model_impedance(double frequency, const Params *p)
{
    double Rs = p->Rs;
    double Rct = p->Rct;
    double Cdl = p->Cdl;

    double w = 2 * PI * frequency;
    double omega_tau = Rct * Cdl * w;
    double denom = 1.0 + omega_tau * omega_tau;

    // Real part: Z_real = Rs + Rct / (1 + (Rct * Cdl * ω)²)
    double z_real = Rs + Rct / denom;

    // Imaginary part: Z_imag = -Rct² * Cdl * ω / (1 + (Rct * Cdl * ω)²)
    double z_imag = -Rct * Rct * Cdl * w / denom;

    return (Complex){z_real, z_imag};
}
void model_derivatives(double frequency, const Params *p, Complex *prime)
{
    double Rs = p->Rs;
    double Rct = p->Rct;
    double Cdl = p->Cdl;

    double w = 2 * PI * frequency;
    double omega_tau = Rct * Cdl * w;
    double omega_tau_sq = omega_tau * omega_tau;
    double denom = 1.0 + omega_tau_sq;
    double denom_sq = denom * denom;

    // ∂Z/∂Rs = 1 + 0j
    prime[0] = (Complex){1.0, 0.0};

    // ∂Z/∂Rct = (1 - (w Cdl)^2 Rct^2) / (1 + (w Rct Cdl)^2)^2 - j (2w Rct Cdl)/(1 + (w Rct Cdl)^2)^2
    double cdl_w_sq = Cdl * Cdl * w * w;
    prime[1] = (Complex){
        (1.0 - cdl_w_sq * Rct * Rct) / denom_sq,
        -2.0 * Rct * Cdl * w / denom_sq};

    // ∂Z/∂Cdl = -2 Rct Cdl w^2 / (1 + (w Rct Cdl)^2)^2 - j (Rct^2 ω) / (1 + (w Rct Cdl)^2)^2
    prime[2] = (Complex){
        -2.0 * Rct * Cdl * w * w / denom_sq,
        -Rct * Rct * w / denom_sq};
}

/* --- Data generation --- */
static void generate_randles_simulation_data(void)
{
    /* True parameter values */
    const double Rs_true = 10.0;   /* 10 ohm */
    const double Rct_true = 100.0; /* 100 ohm */
    const double Cdl_true = 1e-6;  /* 1 uF */

    Params true_params = {
        .Rs = Rs_true,
        .Rct = Rct_true,
        .Cdl = Cdl_true};

    /* Frequency range: 0.01 Hz to 100 kHz (logarithmic) */
    const double freq_min = 0.01;     /* 0.01 Hz */
    const double freq_max = 100000.0; /* 100 kHz */

    /* Generate logarithmic frequency range */
    double log_freq_min = log10(freq_min);
    double log_freq_max = log10(freq_max);

    for (int i = 0; i < N; i++)
    {
        /* Logarithmic frequency spacing */
        double log_freq = log_freq_min + (log_freq_max - log_freq_min) * i / (N - 1);
        frequencies[i] = pow(10.0, log_freq);

        /* Calculate ideal impedance */
        impedances[i] = model_impedance(frequencies[i], &true_params);
    }
    printf("\nData generation complete.\n\n");
}

/* --- Weighting system --- */
void set_weighting_method(WeightingMethod method)
{
    current_weighting = method;

    const char *method_names[] = {
        "Unity (no weighting)",
        "Modulus measured (1/|Z_meas|)",
        "Proportional measured (1/|Re_meas|, 1/|Im_meas|)",
        "Modulus fitted (1/|Z_fit|)",
        "Proportional fitted (1/|Re_fit|, 1/|Im_fit|)",
        "Modulus auto (measured->fitted)",
        "Proportional auto (measured->fitted)"};

    printf("Weighting method set to: %s\n", method_names[method]);
}
static void calculate_weights(Complex z_measured, Complex z_fitted,
                              int use_measured_for_auto,
                              double *weight_real, double *weight_imag)
{
    switch (current_weighting)
    {
    case WEIGHT_UNITY:
        *weight_real = 1.0;
        *weight_imag = 1.0;
        break;

    case WEIGHT_MODULUS_MEAS:
    {
        /* Weight by 1/|Z_measured| */
        double magnitude = sqrt(z_measured.real * z_measured.real +
                                z_measured.imag * z_measured.imag);
        double weight = (magnitude > 1e-12) ? (1.0 / magnitude) : 1e12;
        *weight_real = weight;
        *weight_imag = weight;
        break;
    }

    case WEIGHT_PROPORTIONAL_MEAS:
    {
        /* Weight by 1/|Re_measured| and 1/|Im_measured| separately */
        double real_abs = fabs(z_measured.real);
        double imag_abs = fabs(z_measured.imag);

        *weight_real = (real_abs > 1e-12) ? (1.0 / real_abs) : 1e12;
        *weight_imag = (imag_abs > 1e-12) ? (1.0 / imag_abs) : 1e12;
        break;
    }

    case WEIGHT_MODULUS_FIT:
    {
        /* Weight by 1/|Z_fitted| */
        double magnitude = sqrt(z_fitted.real * z_fitted.real +
                                z_fitted.imag * z_fitted.imag);
        double weight = (magnitude > 1e-12) ? (1.0 / magnitude) : 1e12;
        *weight_real = weight;
        *weight_imag = weight;
        break;
    }

    case WEIGHT_PROPORTIONAL_FIT:
    {
        /* Weight by 1/|Re_fitted| and 1/|Im_fitted| separately */
        double real_abs = fabs(z_fitted.real);
        double imag_abs = fabs(z_fitted.imag);

        *weight_real = (real_abs > 1e-12) ? (1.0 / real_abs) : 1e12;
        *weight_imag = (imag_abs > 1e-12) ? (1.0 / imag_abs) : 1e12;
        break;
    }

    default:
        *weight_real = 1.0;
        *weight_imag = 1.0;
        break;
    }
}

/* --- Optimization functions --- */
void evaluate_residuals(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
    (void)info;
    (void)data;
    (void)m_dat;

    // Simple 3-parameter Randles model
    Params p = {
        .Rs = par[0],
        .Rct = par[1],
        .Cdl = par[2]};

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        /* Calculate separate weighting factors for real and imaginary parts */
        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, current_weighting,
                          &weight_real, &weight_imag);

        /* Apply weighted residuals */
        fvec[i] = weight_real * (z_fitted.real - z_measured.real);     /* Real part */
        fvec[i + N] = weight_imag * (z_fitted.imag - z_measured.imag); /* Imaginary part */
    }
}
void evaluate_jacobian(const double *par, int m_dat, const void *data,
                       double *fjac, int *info)
{
    (void)info;
    (void)data;
    (void)m_dat;

    // Simple 3-parameter Randles model
    Params p = {
        .Rs = par[0],
        .Rct = par[1],
        .Cdl = par[2]};

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        /* Calculate model derivatives */
        Complex derivatives[3];
        model_derivatives(frequencies[i], &p, derivatives);

        /* Calculate separate weighting factors for real and imaginary parts */
        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, current_weighting,
                          &weight_real, &weight_imag);

        /* Jacobian for real part residual (row i) */
        fjac[0 * m_dat + i] = weight_real * derivatives[0].real; /* ∂residual_real/∂Rs */
        fjac[1 * m_dat + i] = weight_real * derivatives[1].real; /* ∂residual_real/∂Rct */
        fjac[2 * m_dat + i] = weight_real * derivatives[2].real; /* ∂residual_real/∂Cdl */

        /* Jacobian for imaginary part residual (row i+N) */
        fjac[0 * m_dat + i + N] = weight_imag * derivatives[0].imag; /* ∂residual_imag/∂Rs */
        fjac[1 * m_dat + i + N] = weight_imag * derivatives[1].imag; /* ∂residual_imag/∂Rct */
        fjac[2 * m_dat + i + N] = weight_imag * derivatives[2].imag; /* ∂residual_imag/∂Cdl */
    }
}
double calculate_chi_squared(const double *par)
{
    Params p = {par[0], par[1], par[2]};
    double chi_sq = 0.0;

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, 0, &weight_real, &weight_imag);

        double residual_real = z_fitted.real - z_measured.real;
        double residual_imag = z_fitted.imag - z_measured.imag;

        chi_sq += weight_real * weight_real * residual_real * residual_real +
                  weight_imag * weight_imag * residual_imag * residual_imag;
    }

    return chi_sq;
}

/* --- Bounds management --- */
lm_bounds_struct *setup_bounds(int n_par, lm_bound_type bound_type, int use_user_scales)
{
    if (n_par <= 0)
    {
        printf("Invalid n_par: %d\n", n_par);
        return NULL;
    }

    lm_bounds_struct *bounds = malloc(sizeof(lm_bounds_struct));
    if (!bounds)
    {
        printf("malloc failed for lm_bounds_struct\n");
        return NULL;
    }

    /* Initialize all pointers to NULL first */
    bounds->lower = NULL;
    bounds->upper = NULL;
    bounds->scales = NULL;
    bounds->bound_type = NULL;

    bounds->lower = malloc(n_par * sizeof(double));
    bounds->upper = malloc(n_par * sizeof(double));
    bounds->bound_type = malloc(n_par * sizeof(int));

    if (use_user_scales)
    {
        bounds->scales = malloc(n_par * sizeof(double));
    }
    else
    {
        bounds->scales = NULL; /* Auto-scaling will be used */
    }

    if (!bounds->lower || !bounds->upper || !bounds->bound_type ||
        (use_user_scales && !bounds->scales))
    {
        printf("Failed to allocate bounds arrays\n");
        free(bounds->lower);
        free(bounds->upper);
        free(bounds->bound_type);
        free(bounds->scales);
        free(bounds);
        return NULL;
    }

    /* Set bounds based on type */
    for (int i = 0; i < n_par; i++)
    {
        switch (bound_type)
        {
        case LM_BOUND_NONE:
            bounds->lower[i] = -1e+6; /* No used */
            bounds->upper[i] = 1e+6;  /* No used */
            bounds->bound_type[i] = LM_BOUND_NONE;
            break;

        case LM_BOUND_LOWER:
            bounds->lower[i] = 1e-12;
            bounds->upper[i] = 1e+6; /* No used */
            bounds->bound_type[i] = LM_BOUND_LOWER;
            break;

        case LM_BOUND_UPPER:
            bounds->lower[i] = 0.0; /* No used */
            bounds->upper[i] = 1e+6;
            bounds->bound_type[i] = LM_BOUND_UPPER;
            break;

        case LM_BOUND_BOTH:
            bounds->lower[i] = 1e-8;
            bounds->upper[i] = 1e+8;
            bounds->bound_type[i] = LM_BOUND_BOTH;
            break;

        case LM_BOUND_FIXED:
            bounds->lower[i] = 0.0; /* Will be set to current parameter value */
            bounds->upper[i] = 0.0; /* Will be set to current parameter value */
            bounds->bound_type[i] = LM_BOUND_FIXED;
            break;

        case LM_BOUND_LOG:
            bounds->lower[i] = 1e-8; /* No used */
            bounds->upper[i] = 1e+8; /* No used */
            bounds->bound_type[i] = LM_BOUND_LOG;
            break;

        default:
            printf("Invalid bound_type: %d\n", bound_type);
            free(bounds->lower);
            free(bounds->upper);
            free(bounds->bound_type);
            free(bounds->scales);
            free(bounds);
            return NULL;
        }
    }

    /* Set user-defined scales if requested */
    if (use_user_scales)
    {
        bounds->scales[0] = 10.0;  /* Rs - resistance scale */
        bounds->scales[1] = 100.0; /* Rct - resistance scale */
        bounds->scales[2] = 1e-6;  /* Cdl - capacitance scale */
    }

    return bounds;
}
void cleanup_bounds(lm_bounds_struct *bounds)
{
    if (bounds)
    {
        free(bounds->lower);
        free(bounds->upper);
        free(bounds->scales);
        free(bounds->bound_type);
        free(bounds);
    }
}

/* --- Testing functions --- */
void test_all_configurations(const double *initial_par)
{
    const int n_par = 3;
    const char *param_names[] = {"Rs", "Rct", "Cdl"};

    /* Test configurations */
    struct test_config
    {
        lm_bound_type bound_type;
        int use_user_scales;
        int use_auto_scales;
        double step_bound;
        const char *config_name;
    } configs[] = {
        {LM_BOUND_NONE, 0, 0, 100.0, "NONE + no scaling + step=100"},
        {LM_BOUND_NONE, 0, 1, 100.0, "NONE + auto scaling + step=100"},
        {LM_BOUND_NONE, 1, 0, 100.0, "NONE + user scaling + step=100"},
        {LM_BOUND_LOWER, 0, 0, 0.001, "LOWER + no scaling + step=0.001"},
        {LM_BOUND_LOWER, 0, 1, 0.001, "LOWER + auto scaling + step=0.001"},
        {LM_BOUND_LOWER, 1, 0, 0.001, "LOWER + user scaling + step=0.001"},
        {LM_BOUND_BOTH, 0, 0, 0.001, "BOTH + no scaling + step=0.001"},
        {LM_BOUND_LOG, 0, 0, 0.001, "LOG + no scaling + step=0.001"},
        {LM_BOUND_LOG, 0, 1, 0.001, "LOG + auto scaling + step=0.001"},
        {LM_BOUND_LOG, 1, 0, 0.001, "LOG + user scaling + step=0.001"}};

    WeightingMethod methods[] = {
        WEIGHT_UNITY,
        WEIGHT_MODULUS_MEAS,
        WEIGHT_PROPORTIONAL_MEAS,
        WEIGHT_MODULUS_FIT,
        WEIGHT_PROPORTIONAL_FIT};

    const char *method_names[] = {
        "Unity",
        "Modulus Measured",
        "Proportional Measured",
        "Modulus Fitted",
        "Proportional Fitted"};

    int num_configs = sizeof(configs) / sizeof(configs[0]);
    int num_methods = sizeof(methods) / sizeof(methods[0]);

    /* Test both Jacobian types */
    const char *jacobian_types[] = {"Numerical", "Analytical"};

    for (int jac_type = 0; jac_type < 2; jac_type++)
    {
        printf("\n=== Testing All Configurations with %s Jacobian ===\n", jacobian_types[jac_type]);
        printf("True parameters:    Rs=10.0, Rct=100.0, Cdl=1e-06\n");
        printf("Initial parameters: Rs=1,    Rct=1,     Cdl=1\n\n");

        /* Results summary table */
        printf("Configuration                          | Weighting Method      | chi2_reduced  | Rs         | Rct        | Cdl        | Outcome\n");
        printf("---------------------------------------|-----------------------|---------------|------------|------------|------------|------------------\n");

        for (int c = 0; c < num_configs; c++)
        {
            /* Setup bounds for current configuration */
            lm_bounds_struct *bounds = setup_bounds(n_par, configs[c].bound_type, configs[c].use_user_scales);
            if (!bounds)
            {
                printf("Error: Failed to allocate bounds for config %d\n", c);
                continue;
            }

            for (int m = 0; m < num_methods; m++)
            {
                /* Reset parameters to initial guess */
                double par[3] = {initial_par[0], initial_par[1], initial_par[2]};
                double par_errors[3] = {0};
                double covar[9] = {0}; /* 3x3 covariance matrix */

                /* Set weighting method */
                current_weighting = methods[m];

                /* Setup control structure */
                lm_control_struct control = lm_control_double;
                control.stepbound = configs[c].step_bound;
                control.patience = 1000;
                control.scale_diag = configs[c].use_auto_scales;
                control.bounds = bounds;
                control.verbosity = 0; // Silent for summary table

                /* Set Jacobian type */
                if (jac_type == 0)
                {
                    control.jacobian = NULL; /* Numerical Jacobian */
                }
                else
                {
                    control.jacobian = evaluate_jacobian; /* Analytical Jacobian */
                }

                control.ftol = 1e-12;
                control.xtol = 1e-12;
                control.gtol = 1e-12;

                /* Run optimization */
                lm_status_struct status;
                lmmin2(n_par, par, par_errors, covar, 2 * N, NULL, NULL, evaluate_residuals, &control, &status);

                /* Calculate degrees of freedom and reduced chi-squared */
                int fixed_params = 0;
                for (int i = 0; i < n_par; i++)
                {
                    if (!bounds->bound_type || bounds->bound_type[i] != LM_BOUND_FIXED)
                    {
                        fixed_params++;
                    }
                }
                int dof = 2 * N - fixed_params;
                double chi2_reduced = (dof > 0) ? (status.fnorm * status.fnorm) / dof : 0.0;

                /* Print summary row */
                printf("%-38s | %-21s | %13.6e | %10.4g | %10.4g | %10.4g | %s\n",
                       configs[c].config_name,
                       method_names[m],
                       chi2_reduced,
                       par[0], par[1], par[2],
                       lm_shortmsg[status.outcome]);

                /* Print standard errors row */
                printf("%-38s | %-21s | %-13s | %10.4g | %10.4g | %10.4g | %s\n",
                       "",
                       "(std errors)",
                       "",
                       par_errors[0], par_errors[1], par_errors[2],
                       "");
            }

            cleanup_bounds(bounds);
            printf("---------------------------------------|-----------------------|---------------|------------|------------|------------|------------------\n");
        }

        printf("\nLegend:\n");
        printf("  NONE:  No bounds (unbounded optimization)\n");
        printf("  LOWER: Lower bounds only [1e-12, inf)\n");
        printf("  UPPER: Upper bounds only (-inf, 1e+6]\n");
        printf("  BOTH:  Both lower and upper bounds [1e-8, 1e+8]\n");
        printf("  FIXED: Fixed parameters (not optimized)\n");
        printf("  LOG:   Logarithmic transformation (always positive)\n");
        printf("  step:  stepbound parameter for trust region\n");
        printf("  auto scaling: automatic parameter scaling by lmmin\n");
        printf("  user scaling: Rs=10, Rct=100, Cdl=1e-6\n\n");
    }
}

/* --- Main function --- */
int main()
{
    // Initial guess for parameters (close to true values)
    double initial_par[3] = {
        1.0, // Rs
        1.0, // Rct
        1.0  // Cdl
    };

    // Generate test data first
    generate_randles_simulation_data();

    // Test all configurations at once
    test_all_configurations(initial_par);

    return 0;
}