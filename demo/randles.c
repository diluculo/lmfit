#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "lmmin.h"

#define N 70
#define PI 3.14159265358979323846

/* --- Complex number structure --- */
typedef struct
{
    double real;
    double imag;
} Complex;

/* --- Basic complex arithmetic --- */

static inline Complex c_add(Complex a, Complex b)
{
    return (Complex){a.real + b.real, a.imag + b.imag};
}

static inline Complex c_sub(Complex a, Complex b)
{
    return (Complex){a.real - b.real, a.imag - b.imag};
}

static inline Complex c_mul(Complex a, Complex b)
{
    return (Complex){
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real};
}

static inline Complex c_mul_scalar(Complex a, double b)
{
    return (Complex){a.real * b, a.imag * b};
}

static inline Complex c_div(Complex a, Complex b)
{
    double denom = b.real * b.real + b.imag * b.imag;
    return (Complex){
        (a.real * b.real + a.imag * b.imag) / denom,
        (a.imag * b.real - a.real * b.imag) / denom};
}

static inline Complex c_inv(Complex a)
{
    double denom = a.real * a.real + a.imag * a.imag;
    return (Complex){a.real / denom, -a.imag / denom};
}

/* --- Complex impedance compositions --- */

static inline Complex z_ser(Complex z1, Complex z2)
{
    return c_add(z1, z2);
}

static inline Complex z_ser_array(int count, const Complex *z)
{
    Complex sum = {0.0, 0.0};
    for (int i = 0; i < count; ++i)
        sum = z_ser(sum, z[i]);
    return sum;
}

static inline Complex z_par(Complex z1, Complex z2)
{
    return c_inv(c_add(c_inv(z1), c_inv(z2)));
}

/* --- Component models --- */

// Resistor: R → Z = R
static inline Complex resistor(double R)
{
    return (Complex){R, 0.0};
}

// Capacitor: C → Z = -j / (ωC)
static inline Complex capacitor(double C, double omega)
{
    return (Complex){0.0, -1.0 / (omega * C)};
}

/* --- Simple Randles circuit parameters --- */
typedef struct
{
    double Rs, Rct, Cdl;
} Params;

/* --- Frequency and measured data --- */

static double frequencies[N];
static Complex impedances[N];

Complex model_impedance(double frequency, const Params *p)
{
    double w = 2 * PI * frequency;

    Complex Z_Rs = resistor(p->Rs);
    Complex Z_Rct = resistor(p->Rct);
    Complex Z_Cdl = capacitor(p->Cdl, w);

    /* Rs in series with (Rct || Cdl) */
    Complex Z_rc = z_par(Z_Rct, Z_Cdl);
    return z_ser(Z_Rs, Z_rc);
}

/* --- Generate ideal Randles circuit data --- */
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

    // Seed random number generator for reproducible noise
    srand(12345);

    for (int i = 0; i < N; i++)
    {
        /* Logarithmic frequency spacing */
        double log_freq = log_freq_min + (log_freq_max - log_freq_min) * i / (N - 1);
        frequencies[i] = pow(10.0, log_freq);

        /* Calculate ideal impedance */
        impedances[i] = model_impedance(frequencies[i], &true_params);

        // /* Add small amount of noise (0.1% of magnitude) */
        // double noise_level = 0.001; /* 0.1% noise */
        // double mag = sqrt(impedances[i].real * impedances[i].real +
        //                   impedances[i].imag * impedances[i].imag);

        // /* Generate random noise */
        // double noise_real = ((double)rand() / RAND_MAX - 0.5) * 2 * noise_level * mag;
        // double noise_imag = ((double)rand() / RAND_MAX - 0.5) * 2 * noise_level * mag;

        // impedances[i].real += noise_real;
        // impedances[i].imag += noise_imag;

        // Print first few and last few points
        // if (i < 3 || i >= N - 3)
        // {
        //     printf("f[%2d] = %8.2e Hz, Z = %8.2f %+8.2fi ohm\n",
        //            i, frequencies[i], impedances[i].real, impedances[i].imag);
        // }
        // else if (i == 3)
        // {
        //     printf("...\n");
        // }
    }
    printf("\nData generation complete.\n\n");
}

/* --- Weighting system --- */

/* Weighting method enumeration */
typedef enum
{
    WEIGHT_UNITY = 0,         /* No weighting (unit weights) */
    WEIGHT_MODULUS_MEAS,      /* Weight by 1/|Z_measured| */
    WEIGHT_PROPORTIONAL_MEAS, /* Weight by 1/|Re|, 1/|Im| separately */
    WEIGHT_MODULUS_FIT,       /* Weight by 1/|Z_fitted| */
    WEIGHT_PROPORTIONAL_FIT,  /* Weight by 1/|Re_fit|, 1/|Im_fit| separately */
    WEIGHT_MODULUS_AUTO,      /* Start with measured, then use fitted */
    WEIGHT_PROPORTIONAL_AUTO  /* Start with measured proportional, then use fitted */
} WeightingMethod;

/* Global weighting method and state */
static WeightingMethod current_weighting = WEIGHT_UNITY; /* Default to no weighting */
static int evaluation_count = 0;                         /* Track evaluations for auto modes */

/* Function to set weighting method */
void set_weighting_method(WeightingMethod method)
{
    current_weighting = method;
    evaluation_count = 0; /* Reset counter when method changes */

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

/* Calculate weighting factors for real and imaginary parts separately */
static void calculate_weights(Complex z_measured, Complex z_fitted,
                              int use_measured_for_auto,
                              double *weight_real, double *weight_imag)
{
    /* For auto modes, decide which data to use */
    WeightingMethod effective_method = current_weighting;
    if (current_weighting == WEIGHT_MODULUS_AUTO)
    {
        effective_method = use_measured_for_auto ? WEIGHT_MODULUS_MEAS : WEIGHT_MODULUS_FIT;
    }
    else if (current_weighting == WEIGHT_PROPORTIONAL_AUTO)
    {
        effective_method = use_measured_for_auto ? WEIGHT_PROPORTIONAL_MEAS : WEIGHT_PROPORTIONAL_FIT;
    }

    switch (effective_method)
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

/* --- Residual function for lmmin --- */

void evaluate_residuals(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
    (void)info;
    (void)data;
    (void)m_dat;

    evaluation_count++;

    /* For auto modes: use measured data for first few evaluations, then switch to fitted */
    int use_measured_for_auto = (evaluation_count <= 3) &&
                                (current_weighting == WEIGHT_MODULUS_AUTO ||
                                 current_weighting == WEIGHT_PROPORTIONAL_AUTO);

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
        calculate_weights(z_measured, z_fitted, use_measured_for_auto,
                          &weight_real, &weight_imag);

        /* Apply weighted residuals */
        fvec[i] = weight_real * (z_fitted.real - z_measured.real);     /* Real part */
        fvec[i + N] = weight_imag * (z_fitted.imag - z_measured.imag); /* Imaginary part */
    }
}

/* Function to calculate chi-squared with current weighting */
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

/* --- Setup bounds for 3 parameters (Rs, Rct, Cdl) --- */
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

/* Function to print detailed frequency-by-frequency comparison */
void print_frequency_comparison(const double *par)
{
    Params p = {par[0], par[1], par[2]};

    printf("\nFrequency comparison (Magnitude & Phase):\n");
    printf(" Freq [Hz]   |Z_meas|      Phase_meas   |Z_fit|       Phase_fit    |Z|_Error%%  Phase_Error\n");
    printf("-----------  ------------  -----------  ------------  -----------  ----------  -----------\n");

    double total_mag_error = 0.0;
    double total_phase_error = 0.0;

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        /* Calculate magnitude and phase */
        double mag_meas = sqrt(z_measured.real * z_measured.real + z_measured.imag * z_measured.imag);
        double phase_meas = atan2(z_measured.imag, z_measured.real) * 180.0 / PI;

        double mag_fit = sqrt(z_fitted.real * z_fitted.real + z_fitted.imag * z_fitted.imag);
        double phase_fit = atan2(z_fitted.imag, z_fitted.real) * 180.0 / PI;

        /* Calculate errors */
        double mag_error_pct = 100.0 * (mag_fit - mag_meas) / mag_meas;
        double phase_error = phase_fit - phase_meas;

        /* Handle phase wraparound */
        if (phase_error > 180.0)
            phase_error -= 360.0;
        if (phase_error < -180.0)
            phase_error += 360.0;

        printf("%11.3f  %9.6e  %11.2f  %9.6e  %11.2f  %10.2f  %11.2f\n",
               frequencies[i], mag_meas, phase_meas, mag_fit, phase_fit,
               mag_error_pct, phase_error);

        /* Accumulate absolute errors for summary */
        total_mag_error += fabs(mag_error_pct);
        total_phase_error += fabs(phase_error);
    }

    printf("-----------  ------------  -----------  ------------  -----------  ----------  -----------\n");
    printf("Average absolute error:                                            %10.2f  %11.2f\n",
           total_mag_error / N, total_phase_error / N);
}

/* Function to test different weighting methods */
void test_weighting_methods(const double *initial_par, lm_bounds_struct *bounds, double step_bound, int use_auto_scales)
{
    const int n_par = 3;
    const char *param_names[] = {"Rs", "Rct", "Cdl"};

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

    int num_methods = sizeof(methods) / sizeof(methods[0]);

    printf("\n=== Testing Different Weighting Methods ===\n");

    for (int m = 0; m < num_methods; m++)
    {
        printf("\n--- %i. Testing %s Weighting ---\n", m + 1, method_names[m]);

        // Reset parameters to initial guess
        double par[3] = {initial_par[0], initial_par[1], initial_par[2]};
        evaluation_count = 0;

        set_weighting_method(methods[m]);

        printf("Initial parameters:\n");
        for (int i = 0; i < n_par; ++i)
            printf("  %-6s = %.10g\n", param_names[i], initial_par[i]);

        double chi_sq = calculate_chi_squared(par);
        printf("initial fnorm: %.6e\n\n", sqrt(chi_sq));

        /* Setup control structure */
        lm_control_struct control = lm_control_double;
        control.stepbound = step_bound;
        control.patience = 1000;
        control.scale_diag = use_auto_scales;
        control.bounds = bounds;
        control.verbosity = 1; // Reduced verbosity for comparison

        control.ftol = 1e-12;
        control.xtol = 1e-12;
        control.gtol = 1e-12;

        lm_status_struct status;
        lmmin2(n_par, par, NULL, NULL, 2 * N, NULL, NULL, evaluate_residuals, &control, &status);
    }
}

/* Function to test different bound types and scaling combinations */
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

    printf("\n=== Testing All Configurations ===\n");
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

            evaluation_count = 0;

            /* Set weighting method */
            current_weighting = methods[m];

            /* Setup control structure */
            lm_control_struct control = lm_control_double;
            control.stepbound = configs[c].step_bound;
            control.patience = 1000;
            control.scale_diag = configs[c].use_auto_scales;
            control.bounds = bounds;
            control.verbosity = 0; // Silent for summary table

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
    printf("  NONE: No bounds (unbounded optimization)\n");
    printf("  BOTH: Both lower and upper bounds [1e-8, 1e+8]\n");
    printf("  LOG:  Logarithmic transformation (always positive)\n");
    printf("  step: stepbound parameter for trust region\n\n");
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