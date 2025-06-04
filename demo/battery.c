#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lmmin.h"

#define N 29
#define PI 3.14159265358979323846

/* --- Complex number structure --- */
typedef struct
{
    double real;
    double imag;
} Complex;

/* --- EIS circuit parameters --- */
typedef struct
{
    double Rstray, L1, Rs, Cf, Rf, Cdl, Rct, Qy1, Qa1;
} Params;

/* --- Global data arrays --- */
static const double frequencies[N] = {
    4000, 3162.28, 2154.43, 1467.8, 1000, 681.292, 464.159, 316.228,
    215.443, 146.78, 100, 68.129, 46.416, 31.623, 21.544, 14.678,
    10, 6.813, 4.641, 3.162, 2.154, 1.468, 1, 0.681,
    0.464, 0.316, 0.215, 0.147, 0.1};
static const Complex impedances[N] = {
    {0.000663089074, 0.00252012886}, {0.000593130962, 0.002012695201}, {0.000514193238, 0.001402412198}, {0.000451633546, 0.000951764771}, {0.000425432787, 0.00062922476}, {0.000420608933, 0.000410356128}, {0.000427075089, 0.000249392344}, {0.000443342498, 0.000133066927}, {0.000466690484, 0.000052738781}, {0.000496198248, 0.000000873319}, {0.000523345566, -0.000028826491}, {0.000551174431, -0.000043343169}, {0.000580768587, -0.000054246608}, {0.000601011623, -0.000057778532}, {0.000616794061, -0.000058350616}, {0.000628636927, -0.000055934551}, {0.000641708479, -0.000056641414}, {0.000652100449, -0.000061967468}, {0.000662929832, -0.000066026}, {0.00067382404, -0.000075871667}, {0.000689054979, -0.0000893174}, {0.000702894694, -0.000108925309}, {0.000719473295, -0.000131383953}, {0.00073985794, -0.000171072749}, {0.000766850441, -0.00021927308}, {0.000801383826, -0.00028480345}, {0.000845546901, -0.000375287209}, {0.000914403545, -0.000494018774}, {0.001018981637, -0.000651889399}};

/* --- Weighting method enumeration --- */
typedef enum
{
    WEIGHT_UNITY = 0,         /* No weighting (unit weights) */
    WEIGHT_MODULUS_MEAS,      /* Weight by 1/|Z_measured| */
    WEIGHT_PROPORTIONAL_MEAS, /* Weight by 1/|Re|, 1/|Im| separately */
    WEIGHT_MODULUS_FIT,       /* Weight by 1/|Z_fitted| */
    WEIGHT_PROPORTIONAL_FIT,  /* Weight by 1/|Re_fit|, 1/|Im_fit| separately */
} WeightingMethod;

/* --- Global weighting state --- */
static WeightingMethod current_weighting = WEIGHT_UNITY;

/* --- Complex arithmetic functions --- */
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

/* --- Impedance composition functions --- */
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

/* --- Component model functions --- */
// Resistor: R → Z = R
static inline Complex resistor(double R)
{
    return (Complex){R, 0.0};
}
// Inductor: L → Z = jωL
static inline Complex inductor(double L, double omega)
{
    return (Complex){0.0, omega * L};
}
// Capacitor: C → Z = -j / (ωC)
static inline Complex capacitor(double C, double omega)
{
    return (Complex){0.0, -1.0 / (omega * C)};
}
// CPE: Q → Z = 1 / (Q * (jω)^a)
static inline Complex constantPhaseElement(double Q, double alpha, double omega)
{
    double magnitude = pow(omega, alpha);
    double phase = alpha * PI / 2.0;

    Complex jw_alpha = {magnitude * cos(phase), magnitude * sin(phase)};
    Complex qz = c_mul_scalar(jw_alpha, Q);
    return c_inv(qz);
}

/* --- EIS model functions --- */
Complex model_impedance(double freq, const Params *p)
{
    double w = 2 * PI * freq;

    Complex Rstray = resistor(p->Rstray);
    Complex L1 = inductor(p->L1, w);
    Complex Rs = resistor(p->Rs);
    Complex Cf = capacitor(p->Cf, w);
    Complex Rf = resistor(p->Rf);
    Complex Cdl = capacitor(p->Cdl, w);
    Complex Rct = resistor(p->Rct);
    Complex Q1 = constantPhaseElement(p->Qy1, p->Qa1, w);

    // Rstray|L1-Rs-Cf|Rf-Cdl|(Rct-Q1)
    Complex Z1 = z_par(Rstray, L1);
    Complex Z2 = z_par(Cf, Rf);
    Complex Z3 = z_par(Cdl, z_ser(Rct, Q1));

    Complex elements[4] = {Z1, Rs, Z2, Z3};
    return z_ser_array(4, elements);
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

    // 9-parameter EIS model
    Params p = {
        .Rstray = par[0],
        .L1 = par[1],
        .Rs = par[2],
        .Cf = par[3],
        .Rf = par[4],
        .Cdl = par[5],
        .Rct = par[6],
        .Qy1 = par[7],
        .Qa1 = par[8]};

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        /* Calculate separate weighting factors for real and imaginary parts */
        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, &weight_real, &weight_imag);

        /* Apply weighted residuals */
        fvec[i] = weight_real * (z_fitted.real - z_measured.real);     /* Real part */
        fvec[i + N] = weight_imag * (z_fitted.imag - z_measured.imag); /* Imaginary part */
    }
}
double calculate_chi_squared(const double *par)
{
    Params p = {par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8]};
    double chi_sq = 0.0;

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, &weight_real, &weight_imag);

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
    if (n_par != 9)
    {
        printf("Invalid n_par: %d (expected 9)\n", n_par);
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
            bounds->lower[i] = -1e+6; /* Not used */
            bounds->upper[i] = 1e+6;  /* Not used */
            bounds->bound_type[i] = LM_BOUND_NONE;
            break;

        case LM_BOUND_LOWER:
            bounds->lower[i] = 1e-12;
            bounds->upper[i] = 1e+6; /* Not used */
            bounds->bound_type[i] = LM_BOUND_LOWER;
            break;

        case LM_BOUND_UPPER:
            bounds->lower[i] = 0.0; /* Not used */
            bounds->upper[i] = 1e+6;
            bounds->bound_type[i] = LM_BOUND_UPPER;
            break;

        case LM_BOUND_BOTH:
            if (i == 8) /* Qa1 (CPE exponent): must be between 0 and 1 */
            {
                bounds->lower[i] = 1e-8;
                bounds->upper[i] = 1.0 - 1e-8;
            }
            else
            {
                bounds->lower[i] = 1e-8;
                bounds->upper[i] = 1e+8;
            }
            bounds->bound_type[i] = LM_BOUND_BOTH;
            break;

        case LM_BOUND_FIXED:
            bounds->lower[i] = 0.0; /* Will be set to current parameter value */
            bounds->upper[i] = 0.0; /* Will be set to current parameter value */
            bounds->bound_type[i] = LM_BOUND_FIXED;
            break;

        case LM_BOUND_LOG:
            bounds->lower[i] = 1e-8; /* Not used */
            bounds->upper[i] = 1e+8; /* Not used */
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
        bounds->scales[0] = 1e-3;   /* Rstray - resistance scale */
        bounds->scales[1] = 1e-6;   /* L1 - inductance scale */
        bounds->scales[2] = 1e-3;   /* Rs - resistance scale */
        bounds->scales[3] = 1.0;    /* Cf - capacitance scale */
        bounds->scales[4] = 1e-3;   /* Rf - resistance scale */
        bounds->scales[5] = 1.0;    /* Cdl - capacitance scale */
        bounds->scales[6] = 1e-3;   /* Rct - resistance scale */
        bounds->scales[7] = 1000.0; /* Qy1 - CPE magnitude scale */
        bounds->scales[8] = 0.5;    /* Qa1 - CPE exponent scale (0-1) */
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
    const int n_par = 9;
    const char *param_names[] = {"Rstray", "L1", "Rs", "Cf", "Rf", "Cdl", "Rct", "Qy1", "Qa1"};

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
    printf("Initial parameters: Rstray=0.0001, L1=2e-7, Rs=0.00055, Cf=2.5, Rf=0.00028, Cdl=18.0, Rct=0.00022, Qy1=6000, Qa1=0.7\n\n");

    /* Results summary table */
    printf("Configuration                          | Weighting Method      | chi2_reduced  | Rstray     | L1         | Rs         | Cf         | Rf         | Cdl        | Rct        | Qy1        | Qa1        | Status\n");
    printf("---------------------------------------|-----------------------|---------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------------\n");

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
            double par[9] = {initial_par[0], initial_par[1], initial_par[2], initial_par[3], initial_par[4],
                             initial_par[5], initial_par[6], initial_par[7], initial_par[8]};
            double par_errors[9] = {0};
            double covar[81] = {0}; /* 9x9 covariance matrix */

            /* Set weighting method */
            current_weighting = methods[m];

            /* Setup control structure */
            lm_control_struct control = lm_control_double;
            control.stepbound = configs[c].step_bound;
            control.patience = 1000;
            control.scale_diag = configs[c].use_auto_scales;
            control.bounds = bounds;
            control.verbosity = 0; // Silent for summary table

            control.ftol = 1e-8;
            control.xtol = 1e-8;
            control.gtol = 1e-8;

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

            /* Print summary row with all fitted parameters */
            printf("%-38s | %-21s | %13.6e | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %s\n",
                   configs[c].config_name,
                   method_names[m],
                   chi2_reduced,
                   par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8],
                   lm_shortmsg[status.outcome]);

            /* Print standard errors row */
            printf("%-38s | %-21s | %13s | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %10.4g | %s\n",
                   "",
                   "(std errors)",
                   "",
                   par_errors[0], par_errors[1], par_errors[2], par_errors[3],
                   par_errors[4], par_errors[5], par_errors[6], par_errors[7], par_errors[8],
                   "");
        }

        cleanup_bounds(bounds);
        printf("---------------------------------------|-----------------------|---------------|------------|------------|------------|------------|------------|------------|------------|------------|------------|------------------\n");
    }

    printf("\nLegend:\n");
    printf("  NONE:  No bounds (unbounded optimization)\n");
    printf("  LOWER: Lower bounds only [1e-12, inf)\n");
    printf("  UPPER: Upper bounds only (-inf, 1e+6]\n");
    printf("  BOTH:  Both bounds [1e-8, 1e+8], Qa1: [1e-8, 1-1e-8]\n");
    printf("  FIXED: Fixed parameters (not optimized)\n");
    printf("  LOG:   Logarithmic transformation (always positive)\n");
    printf("  step:  stepbound parameter for trust region\n");
    printf("  scaling: auto=lmmin internal, user=custom per parameter\n\n");
}

/* --- Main function --- */
int main()
{
    // Initial guess for parameters
    double initial_par[9] = {
        0.0001,  // Rstray
        2e-7,    // L1
        0.00055, // Rs
        2.5,     // Cf
        0.00028, // Rf
        18.0,    // Cdl
        0.00022, // Rct
        6000,    // Qy1
        0.7      // Qa1
    };

    // Test all configurations at once
    test_all_configurations(initial_par);

    return 0;
}