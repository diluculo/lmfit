#ifndef EIS_FRAMEWORK_H
#define EIS_FRAMEWORK_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include "lmmin.h"

#define PI 3.14159265358979323846
#define MAX_PARAMS 20
#define MAX_FREQUENCIES 1000

/* --- Complex number structure --- */
typedef struct
{
    double real;
    double imag;
} Complex;

/* --- Data structure for frequency-impedance pairs --- */
typedef struct
{
    int n_points;
    double *frequencies;
    Complex *impedances;
    char *data_source;
} EISData;

/* --- Generic parameter structure --- */
typedef struct
{
    int n_params;
    double *values;
    char **names;
    double *lower_bounds;
    double *upper_bounds;
    int *bound_types;
    double *scales;
} Parameters;

/* --- Model function pointer --- */
typedef Complex (*ModelFunction)(double frequency, const Parameters *params);

/* --- Weighting method enumeration --- */
typedef enum
{
    WEIGHT_UNITY = 0,
    WEIGHT_MODULUS_MEAS,
    WEIGHT_PROPORTIONAL_MEAS,
    WEIGHT_MODULUS_FIT,
    WEIGHT_PROPORTIONAL_FIT,
    WEIGHT_MODULUS_AUTO,
    WEIGHT_PROPORTIONAL_AUTO
} WeightingMethod;

/* --- Fitting context structure --- */
typedef struct
{
    EISData *data;
    ModelFunction model_func;
    Parameters *params;
    WeightingMethod weighting;
    int evaluation_count;
    double *residuals_array;
} FittingContext;

/* --- Function declarations --- */

// Complex arithmetic
Complex c_add(Complex a, Complex b);
Complex c_sub(Complex a, Complex b);
Complex c_mul(Complex a, Complex b);
Complex c_mul_scalar(Complex a, double b);
Complex c_div(Complex a, Complex b);
Complex c_inv(Complex a);

// Impedance compositions
Complex z_ser(Complex z1, Complex z2);
Complex z_ser_array(int count, const Complex *z);
Complex z_par(Complex z1, Complex z2);

// Basic components
Complex resistor(double R);
Complex inductor(double L, double omega);
Complex capacitor(double C, double omega);
Complex constantPhaseElement(double Q, double alpha, double omega);
Complex warburg(double sigma, double omega);

// Data management
EISData *create_eis_data(int n_points);
void free_eis_data(EISData *data);
int load_data_from_file(const char *filename, EISData *data);
int load_data_from_arrays(double *freq, Complex *imp, int n_points, EISData *data);

// Parameter management
Parameters *create_parameters(int n_params);
void free_parameters(Parameters *params);
int set_parameter(Parameters *params, int index, const char *name, double value,
                  double lower, double upper, int bound_type, double scale);

// Fitting context
FittingContext *create_fitting_context(EISData *data, ModelFunction model_func,
                                       Parameters *params);
void free_fitting_context(FittingContext *ctx);

// Weighting functions
void set_weighting_method(FittingContext *ctx, WeightingMethod method);
void calculate_weights(Complex z_measured, Complex z_fitted,
                       WeightingMethod method, int use_measured_for_auto,
                       double *weight_real, double *weight_imag);

// Fitting functions
void evaluate_residuals(const double *par, int m_dat, const void *data,
                        double *fvec, int *info);
int fit_eis_model(FittingContext *ctx, lm_control_struct *control);

// Analysis and output
double calculate_chi_squared(FittingContext *ctx, const double *par);
void print_fit_results(FittingContext *ctx, const double *par,
                       const lm_status_struct *status);
void print_frequency_comparison(FittingContext *ctx, const double *par);
void export_fit_results(FittingContext *ctx, const double *par,
                        const char *filename);

// Bounds management
lm_bounds_struct *setup_bounds_from_parameters(const Parameters *params);
void cleanup_bounds(lm_bounds_struct *bounds);

void print_parameter_sensitivity(FittingContext *ctx, const double *par);

#endif // EIS_FRAMEWORK_H