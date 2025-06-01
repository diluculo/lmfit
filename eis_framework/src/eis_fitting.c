#include "eis_framework.h"

/* --- Global fitting context (for lmmin callback) --- */
static FittingContext *global_ctx = NULL;

/* --- Fitting context management --- */
FittingContext *create_fitting_context(EISData *data, ModelFunction model_func,
                                       Parameters *params)
{
    if (!data || !model_func || !params)
    {
        printf("Invalid arguments to create_fitting_context\n");
        return NULL;
    }

    FittingContext *ctx = malloc(sizeof(FittingContext));
    if (!ctx)
    {
        printf("Failed to allocate FittingContext\n");
        return NULL;
    }

    ctx->data = data;
    ctx->model_func = model_func;
    ctx->params = params;
    ctx->weighting = WEIGHT_MODULUS_AUTO;
    ctx->evaluation_count = 0;

    // Allocate residuals array
    ctx->residuals_array = malloc(2 * data->n_points * sizeof(double));
    if (!ctx->residuals_array)
    {
        printf("Failed to allocate residuals array\n");
        free(ctx);
        return NULL;
    }

    return ctx;
}

void free_fitting_context(FittingContext *ctx)
{
    if (ctx)
    {
        free(ctx->residuals_array);
        free(ctx);
    }
}

/* --- Weighting functions --- */
void set_weighting_method(FittingContext *ctx, WeightingMethod method)
{
    if (!ctx)
        return;

    ctx->weighting = method;
    ctx->evaluation_count = 0;

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

void calculate_weights(Complex z_measured, Complex z_fitted,
                       WeightingMethod method, int use_measured_for_auto,
                       double *weight_real, double *weight_imag)
{

    WeightingMethod effective_method = method;
    if (method == WEIGHT_MODULUS_AUTO)
    {
        effective_method = use_measured_for_auto ? WEIGHT_MODULUS_MEAS : WEIGHT_MODULUS_FIT;
    }
    else if (method == WEIGHT_PROPORTIONAL_AUTO)
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
        double magnitude_sq = z_measured.real * z_measured.real +
                              z_measured.imag * z_measured.imag;
        double weight = (magnitude_sq > 1e-24) ? (1.0 / magnitude_sq) : 1e24;
        *weight_real = weight;
        *weight_imag = weight;
        break;
    }

    case WEIGHT_PROPORTIONAL_MEAS:
    {
        double real_sq = z_measured.real * z_measured.real;
        double imag_sq = z_measured.imag * z_measured.imag;
        *weight_real = (real_sq > 1e-24) ? (1.0 / real_sq) : 1e24;
        *weight_imag = (imag_sq > 1e-24) ? (1.0 / imag_sq) : 1e24;
        break;
    }

    case WEIGHT_MODULUS_FIT:
    {
        double magnitude_sq = z_fitted.real * z_fitted.real +
                              z_fitted.imag * z_fitted.imag;
        double weight = (magnitude_sq > 1e-24) ? (1.0 / magnitude_sq) : 1e24;
        *weight_real = weight;
        *weight_imag = weight;
        break;
    }

    case WEIGHT_PROPORTIONAL_FIT:
    {
        double real_sq = z_fitted.real * z_fitted.real;
        double imag_sq = z_fitted.imag * z_fitted.imag;
        *weight_real = (real_sq > 1e-24) ? (1.0 / real_sq) : 1e24;
        *weight_imag = (imag_sq > 1e-24) ? (1.0 / imag_sq) : 1e24;
        break;
    }

    default:
        *weight_real = 1.0;
        *weight_imag = 1.0;
        break;
    }
}

/* --- Residual evaluation function --- */
void evaluate_residuals(const double *par, int m_dat, const void *data,
                        double *fvec, int *info)
{
    (void)data; // We use global context instead
    (void)info;

    if (!global_ctx)
    {
        printf("Error: Global fitting context not set\n");
        return;
    }

    FittingContext *ctx = global_ctx;
    ctx->evaluation_count++;

    int use_measured_for_auto = (ctx->evaluation_count == 1) &&
                                (ctx->weighting == WEIGHT_MODULUS_AUTO ||
                                 ctx->weighting == WEIGHT_PROPORTIONAL_AUTO);

    // Update parameter values
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        ctx->params->values[i] = par[i];
    }

    for (int i = 0; i < ctx->data->n_points; i++)
    {
        Complex z_fitted = ctx->model_func(ctx->data->frequencies[i], ctx->params);
        Complex z_measured = ctx->data->impedances[i];

        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, ctx->weighting,
                          use_measured_for_auto, &weight_real, &weight_imag);

        fvec[i] = weight_real * (z_fitted.real - z_measured.real);
        fvec[i + ctx->data->n_points] = weight_imag * (z_fitted.imag - z_measured.imag);
    }
}

/* --- Main fitting function --- */
int fit_eis_model(FittingContext *ctx, lm_control_struct *control)
{
    if (!ctx || !control)
    {
        printf("Invalid arguments to fit_eis_model\n");
        return -1;
    }

    // Set global context for callback
    global_ctx = ctx;

    // Setup bounds
    lm_bounds_struct *bounds = setup_bounds_from_parameters(ctx->params);
    if (!bounds)
    {
        printf("Failed to setup parameter bounds\n");
        return -1;
    }

    control->bounds = bounds;

    // Setup initial values array
    double *par = malloc(ctx->params->n_params * sizeof(double));
    if (!par)
    {
        printf("Failed to allocate parameter array\n");
        cleanup_bounds(bounds);
        return -1;
    }

    for (int i = 0; i < ctx->params->n_params; i++)
    {
        par[i] = ctx->params->values[i];
    }

    // Prepare data array for lmmin
    for (int i = 0; i < ctx->data->n_points; i++)
    {
        ctx->residuals_array[i] = ctx->data->impedances[i].real;
        ctx->residuals_array[i + ctx->data->n_points] = ctx->data->impedances[i].imag;
    }

    lm_status_struct status;

    // Perform fitting
    lmmin2(ctx->params->n_params, par, NULL, NULL,
           2 * ctx->data->n_points, ctx->residuals_array, NULL,
           evaluate_residuals, control, &status);

    // Update parameter values with fitted results
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        ctx->params->values[i] = par[i];
    }

    // Print results
    print_fit_results(ctx, par, &status);

    // Cleanup
    free(par);
    cleanup_bounds(bounds);
    global_ctx = NULL;

    return status.outcome;
}

/* --- Bounds management --- */
lm_bounds_struct *setup_bounds_from_parameters(const Parameters *params)
{
    if (!params)
        return NULL;

    lm_bounds_struct *bounds = malloc(sizeof(lm_bounds_struct));
    if (!bounds)
    {
        printf("malloc failed for lm_bounds_struct\n");
        return NULL;
    }

    bounds->lower = malloc(params->n_params * sizeof(double));
    bounds->upper = malloc(params->n_params * sizeof(double));
    bounds->scales = malloc(params->n_params * sizeof(double));
    bounds->bound_type = malloc(params->n_params * sizeof(int));

    if (!bounds->lower || !bounds->upper || !bounds->scales || !bounds->bound_type)
    {
        printf("Failed to allocate bounds arrays\n");
        cleanup_bounds(bounds);
        return NULL;
    }

    for (int i = 0; i < params->n_params; i++)
    {
        bounds->lower[i] = params->lower_bounds[i];
        bounds->upper[i] = params->upper_bounds[i];
        bounds->scales[i] = params->scales[i];
        bounds->bound_type[i] = params->bound_types[i];
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