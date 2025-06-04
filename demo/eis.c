#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "lmmin.h"

#define N 29
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
    // Z_CPE = 1 / (Q * (jω)^α)
    // (jω)^α = ω^α * e^(j*α*π/2) = ω^α * (cos(α*π/2) + j*sin(α*π/2))

    double magnitude = pow(omega, alpha);
    double phase = alpha * PI / 2.0;

    Complex jw_alpha = {magnitude * cos(phase), magnitude * sin(phase)};

    // Z = 1 / (Q * (jω)^α)
    Complex qz = c_mul_scalar(jw_alpha, Q);
    return c_inv(qz);
}

/* --- Parameters --- */

typedef struct
{
    double Rstray, L1, Rs, Cf, Rf, Cdl, Rct, Qy1, Qa1;
} Params;

/* --- Frequency and measured data --- */

static const double frequencies[N] = {
    4000, 3162.28, 2154.43, 1467.8, 1000, 681.292, 464.159, 316.228,
    215.443, 146.78, 100, 68.129, 46.416, 31.623, 21.544, 14.678,
    10, 6.813, 4.641, 3.162, 2.154, 1.468, 1, 0.681,
    0.464, 0.316, 0.215, 0.147, 0.1};

static const Complex impedances[N] = {
    {0.000663089074, 0.00252012886}, {0.000593130962, 0.002012695201}, {0.000514193238, 0.001402412198}, {0.000451633546, 0.000951764771}, {0.000425432787, 0.00062922476}, {0.000420608933, 0.000410356128}, {0.000427075089, 0.000249392344}, {0.000443342498, 0.000133066927}, {0.000466690484, 0.000052738781}, {0.000496198248, 0.000000873319}, {0.000523345566, -0.000028826491}, {0.000551174431, -0.000043343169}, {0.000580768587, -0.000054246608}, {0.000601011623, -0.000057778532}, {0.000616794061, -0.000058350616}, {0.000628636927, -0.000055934551}, {0.000641708479, -0.000056641414}, {0.000652100449, -0.000061967468}, {0.000662929832, -0.000066026}, {0.00067382404, -0.000075871667}, {0.000689054979, -0.0000893174}, {0.000702894694, -0.000108925309}, {0.000719473295, -0.000131383953}, {0.00073985794, -0.000171072749}, {0.000766850441, -0.00021927308}, {0.000801383826, -0.00028480345}, {0.000845546901, -0.000375287209}, {0.000914403545, -0.000494018774}, {0.001018981637, -0.000651889399}};

/* --- Model impedance function --- */

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

/* --- Residual function for lmmin --- */

static double y[2 * N];

void initialize_y_array()
{
    for (int i = 0; i < N; ++i)
    {
        y[i] = impedances[i].real;
        y[i + N] = impedances[i].imag;
    }
}

/* Weighting method enumeration */
typedef enum
{
    WEIGHT_UNITY = 0,         /* No weighting (unit weights) */
    WEIGHT_MODULUS_MEAS,      /* Weight by 1/|Z_measured|^2 */
    WEIGHT_PROPORTIONAL_MEAS, /* Weight by 1/|Re|^2 and 1/|Im|^2 separately */
    WEIGHT_MODULUS_FIT,       /* Weight by 1/|Z_fitted|^2 */
    WEIGHT_PROPORTIONAL_FIT,  /* Weight by 1/|Re_fit|^2 and 1/|Im_fit|^2 separately */
    WEIGHT_MODULUS_AUTO,      /* Start with measured, then use fitted */
    WEIGHT_PROPORTIONAL_AUTO  /* Start with measured proportional, then use fitted */
} WeightingMethod;

/* Global weighting method and state */
static WeightingMethod current_weighting = WEIGHT_MODULUS_AUTO; /* Default to auto mode */
static int evaluation_count = 1;                                /* Flag to track first evaluation for auto modes */

/* Function to set weighting method */
void set_weighting_method(WeightingMethod method)
{
    current_weighting = method;
    evaluation_count = 1; /* Reset flag when method changes */

    const char *method_names[] = {
        "Unity (no weighting)",
        "Modulus measured (1/|Z_meas|^2)",
        "Proportional measured (1/|Re_meas|^2, 1/|Im_meas|^2)",
        "Modulus fitted (1/|Z_fit|^2)",
        "Proportional fitted (1/|Re_fit|^2, 1/|Im_fit|^2)",
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
        double magnitude_sq = z_measured.real * z_measured.real +
                              z_measured.imag * z_measured.imag;
        double magnitude = sqrt(magnitude_sq);
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

        *weight_real = (real_abs > 1e-12) ? (1.0 / real_abs) : 1e12; // 1/|Re|
        *weight_imag = (imag_abs > 1e-12) ? (1.0 / imag_abs) : 1e12; // 1/|Im|
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

/* Enhanced residual function with separate real/imaginary weighting */
void evaluate_residuals(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
    (void)info;
    (void)data;

    evaluation_count++;

    int use_measured_for_auto = (evaluation_count == 1) &&
                                (current_weighting == WEIGHT_MODULUS_AUTO ||
                                 current_weighting == WEIGHT_PROPORTIONAL_AUTO);

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
        calculate_weights(z_measured, z_fitted, use_measured_for_auto, &weight_real, &weight_imag);

        /* Apply weighted residuals */
        fvec[i] = weight_real * (z_fitted.real - z_measured.real);     /* Real part */
        fvec[i + N] = weight_imag * (z_fitted.imag - z_measured.imag); /* Imaginary part */
    }
}

/* Function to calculate chi-squared with current weighting */
double calculate_chi_squared(const double *par)
{
    Params p = {par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8]};

    double chi_sq = 0.0;

    for (int i = 0; i < N; ++i)
    {
        Complex z_fitted = model_impedance(frequencies[i], &p);
        Complex z_measured = impedances[i];

        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, current_weighting, &weight_real, &weight_imag);

        double residual_real = z_fitted.real - z_measured.real;
        double residual_imag = z_fitted.imag - z_measured.imag;

        chi_sq += weight_real * weight_real * residual_real * residual_real +
                  weight_imag * weight_imag * residual_imag * residual_imag;
    }

    return chi_sq;
}

/* --- Setup bounds for all parameters (0 to infinity) --- */
lm_bounds_struct *setup_bounds(int n_par)
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
        printf("errno = %d\n", errno);
        return NULL;
    }

    /* Initialize all pointers to NULL first */
    bounds->lower = NULL;
    bounds->upper = NULL;
    bounds->scales = NULL;
    bounds->bound_type = NULL;

    bounds->lower = malloc(n_par * sizeof(double));
    if (!bounds->lower)
    {
        printf("Failed to allocate lower array\n");
        free(bounds);
        return NULL;
    }

    bounds->upper = malloc(n_par * sizeof(double));
    if (!bounds->upper)
    {
        printf("Failed to allocate upper array\n");
        free(bounds->lower);
        free(bounds);
        return NULL;
    }

    bounds->scales = malloc(n_par * sizeof(double));
    if (!bounds->scales)
    {
        printf("Failed to allocate scales array\n");
        free(bounds->lower);
        free(bounds->upper);
        free(bounds);
        return NULL;
    }

    bounds->bound_type = malloc(n_par * sizeof(int));
    if (!bounds->bound_type)
    {
        printf("Failed to allocate bound_type array\n");
        free(bounds->lower);
        free(bounds->upper);
        free(bounds->scales);
        free(bounds);
        return NULL;
    }

    /* Set lower bounds to small positive values and appropriate scales */
    for (int i = 0; i < n_par; i++)
    {
        bounds->lower[i] = 0.0;                 /* Small positive lower bound */
        bounds->bound_type[i] = LM_BOUND_LOWER; /* Lower bound only */
    }

    /* Special case for Qa1 (CPE exponent): must be between 0 and 1 */
    // bounds->lower[8] = 0.0;
    // bounds->upper[8] = 1.0;
    // bounds->bound_type[8] = LM_BOUND_BOTH; /* Both bounds */

    /* Set appropriate scales for different parameter types */
    bounds->scales[0] = 1e-3;   /* Rstray - resistance scale */
    bounds->scales[1] = 1e-6;   /* L1 - inductance scale */
    bounds->scales[2] = 1e-3;   /* Rs - resistance scale */
    bounds->scales[3] = 1.0;    /* Cf - capacitance scale */
    bounds->scales[4] = 1e-3;   /* Rf - resistance scale */
    bounds->scales[5] = 1.0;    /* Cdl - capacitance scale */
    bounds->scales[6] = 1e-3;   /* Rct - resistance scale */
    bounds->scales[7] = 1000.0; /* Qy1 - CPE magnitude scale */
    bounds->scales[8] = 0.5;    /* Qa1 - CPE exponent scale (0-1) */

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
    Params p = {par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8]};

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

/* --- Main function --- */

int main()
{
    const int n_par = 9;

    /* Parameter names for better output */
    const char *param_names[] = {
        "Rstray", "L1", "Rs", "Cf", "Rf", "Cdl", "Rct", "Qy1", "Qa1"};

    // Initial guess for parameters
    double par[9] = {
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

    initialize_y_array();

    /* For automatic selection: */
    set_weighting_method(WEIGHT_MODULUS_MEAS);

    /* Setup bounds for all parameters (0 to infinity) */
    lm_bounds_struct *bounds = setup_bounds(n_par);
    if (!bounds)
    {
        fprintf(stderr, "Error: Failed to allocate bounds\n");
        return 1;
    }

    printf("Starting bounded optimization...\n\n");

    /* Setup control structure with bounds */
    lm_control_struct control = lm_control_double;
    control.stepbound = 10.0; /* Initial step bound */
    control.patience = 1000;  /* Maximum number of function evaluations */
    control.scale_diag = 0;   /* Rescale variables internally */
    control.bounds = bounds;  /* Apply bounds */
    control.verbosity = 3;    /* Print more information */

    /* Adjust tolerances for better convergence with bounds */
    control.ftol = 1e-15;
    control.xtol = 1e-15;
    control.gtol = 1e-15;

    lm_status_struct status;

    lmmin2(n_par, par, NULL, NULL, 2 * N, y, NULL, evaluate_residuals, &control, &status);

    printf("\nOptimization completed with status: %s\n", lm_infmsg[status.outcome]);
    printf("Number of function evaluations: %d\n", status.nfev);
    printf("Final residual norm: %.6e\n\n", status.fnorm);

    printf("Fitted parameters:\n");
    for (int i = 0; i < n_par; ++i)
        printf("  %-6s = %.10g\n", param_names[i], par[i]);

    /* Calculate final model fit quality */
    double total_error = 0.0;
    for (int i = 0; i < N; ++i)
    {
        Params p = {par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8]};
        Complex z_model = model_impedance(frequencies[i], &p);
        Complex z_meas = impedances[i];

        double error_real = z_model.real - z_meas.real;
        double error_imag = z_model.imag - z_meas.imag;
        total_error += error_real * error_real + error_imag * error_imag;
    }

    printf("Total squared error: %.6e\n\n", total_error);

    print_frequency_comparison(par);

    cleanup_bounds(bounds);
    return 0;
}
