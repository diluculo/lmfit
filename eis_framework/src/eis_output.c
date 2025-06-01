#include "eis_framework.h"

/* --- Utility function for printing separators --- */
static void print_separator(char c, int count)
{
    for (int i = 0; i < count; i++)
    {
        printf("%c", c);
    }
    printf("\n");
}

/* --- Calculate chi-squared --- */
double calculate_chi_squared(FittingContext *ctx, const double *par)
{
    if (!ctx || !par)
        return -1.0;

    double chi_sq = 0.0;

    // Update parameters
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        ctx->params->values[i] = par[i];
    }

    for (int i = 0; i < ctx->data->n_points; i++)
    {
        Complex z_fitted = ctx->model_func(ctx->data->frequencies[i], ctx->params);
        Complex z_measured = ctx->data->impedances[i];

        double weight_real, weight_imag;
        calculate_weights(z_measured, z_fitted, ctx->weighting, 0,
                          &weight_real, &weight_imag);

        double residual_real = z_fitted.real - z_measured.real;
        double residual_imag = z_fitted.imag - z_measured.imag;

        chi_sq += weight_real * weight_real * residual_real * residual_real +
                  weight_imag * weight_imag * residual_imag * residual_imag;
    }

    return chi_sq;
}

/* --- Print fit results --- */
void print_fit_results(FittingContext *ctx, const double *par,
                       const lm_status_struct *status)
{
    if (!ctx || !par || !status)
        return;

    printf("\n");
    print_separator('=', 60);
    printf("EIS MODEL FITTING RESULTS\n");
    print_separator('=', 60);

    printf("Optimization Status: %s\n", lm_infmsg[status->outcome]);
    printf("Function evaluations: %d\n", status->nfev);
    printf("Final residual norm: %.6e\n", status->fnorm);

    if (ctx->data->data_source)
    {
        printf("Data source: %s\n", ctx->data->data_source);
    }
    printf("Data points: %d\n", ctx->data->n_points);

    printf("\nFitted Parameters:\n");
    print_separator('-', 40);
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        printf("%-10s = %12.6e", ctx->params->names[i], par[i]);

        // Show bounds
        if (ctx->params->bound_types[i] == LM_BOUND_BOTH)
        {
            printf("  [%.2e, %.2e]",
                   ctx->params->lower_bounds[i], ctx->params->upper_bounds[i]);
        }
        else if (ctx->params->bound_types[i] == LM_BOUND_LOWER)
        {
            printf("  [%.2e, +inf)", ctx->params->lower_bounds[i]);
        }
        else if (ctx->params->bound_types[i] == LM_BOUND_UPPER)
        {
            printf("  (-inf, %.2e]", ctx->params->upper_bounds[i]);
        }
        printf("\n");
    }

    double chi_sq = calculate_chi_squared(ctx, par);
    printf("\nGoodness of fit:\n");
    print_separator('-', 20);
    printf("Chi-squared: %.6e\n", chi_sq);
    printf("Reduced chi-squared: %.6e\n",
           chi_sq / (2 * ctx->data->n_points - ctx->params->n_params));

    printf("\n");
}

/* --- Print frequency comparison --- */
void print_frequency_comparison(FittingContext *ctx, const double *par)
{
    if (!ctx || !par)
        return;

    // Update parameters
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        ctx->params->values[i] = par[i];
    }

    printf("\nFrequency-by-frequency comparison:\n");
    print_separator('=', 80);
    printf("%11s %12s %12s %12s %12s %10s %12s\n",
           "Freq [Hz]", "|Z_meas|", "Phase_meas", "|Z_fit|", "Phase_fit",
           "|Z|_Err%", "Phase_Err");
    print_separator('-', 80);

    double total_mag_error = 0.0;
    double total_phase_error = 0.0;
    double max_mag_error = 0.0;
    double max_phase_error = 0.0;

    for (int i = 0; i < ctx->data->n_points; i++)
    {
        Complex z_fitted = ctx->model_func(ctx->data->frequencies[i], ctx->params);
        Complex z_measured = ctx->data->impedances[i];

        // Calculate magnitude and phase
        double mag_meas = sqrt(z_measured.real * z_measured.real +
                               z_measured.imag * z_measured.imag);
        double phase_meas = atan2(z_measured.imag, z_measured.real) * 180.0 / PI;

        double mag_fit = sqrt(z_fitted.real * z_fitted.real +
                              z_fitted.imag * z_fitted.imag);
        double phase_fit = atan2(z_fitted.imag, z_fitted.real) * 180.0 / PI;

        // Calculate errors
        double mag_error_pct = 100.0 * (mag_fit - mag_meas) / mag_meas;
        double phase_error = phase_fit - phase_meas;

        // Handle phase wraparound
        while (phase_error > 180.0)
            phase_error -= 360.0;
        while (phase_error < -180.0)
            phase_error += 360.0;

        printf("%11.3e %12.6e %12.2f %12.6e %12.2f %10.2f %12.2f\n",
               ctx->data->frequencies[i], mag_meas, phase_meas,
               mag_fit, phase_fit, mag_error_pct, phase_error);

        // Accumulate statistics
        total_mag_error += fabs(mag_error_pct);
        total_phase_error += fabs(phase_error);
        if (fabs(mag_error_pct) > max_mag_error)
            max_mag_error = fabs(mag_error_pct);
        if (fabs(phase_error) > max_phase_error)
            max_phase_error = fabs(phase_error);
    }

    print_separator('-', 80);
    printf("Statistics:\n");
    printf("  Average |Z| error:     %10.2f%%\n", total_mag_error / ctx->data->n_points);
    printf("  Average phase error:   %10.2f degrees\n", total_phase_error / ctx->data->n_points);
    printf("  Maximum |Z| error:     %10.2f%%\n", max_mag_error);
    printf("  Maximum phase error:   %10.2f degrees\n", max_phase_error);
    printf("\n");
}

/* --- Export fit results to file --- */
void export_fit_results(FittingContext *ctx, const double *par, const char *filename)
{
    if (!ctx || !par || !filename)
        return;

    FILE *file = fopen(filename, "w");
    if (!file)
    {
        printf("Cannot create output file: %s\n", filename);
        return;
    }

    // Update parameters
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        ctx->params->values[i] = par[i];
    }

    // Write header
    fprintf(file, "# EIS Model Fitting Results\n");
    fprintf(file, "# Generated by EIS Framework\n");
    if (ctx->data->data_source)
    {
        fprintf(file, "# Data source: %s\n", ctx->data->data_source);
    }
    fprintf(file, "# Data points: %d\n", ctx->data->n_points);
    fprintf(file, "#\n");

    // Write parameters
    fprintf(file, "# Fitted Parameters:\n");
    for (int i = 0; i < ctx->params->n_params; i++)
    {
        fprintf(file, "# %s = %.10e\n", ctx->params->names[i], par[i]);
    }
    fprintf(file, "#\n");

    // Write data comparison
    fprintf(file, "# Frequency[Hz] Z_meas_real Z_meas_imag Z_fit_real Z_fit_imag ");
    fprintf(file, "|Z_meas| Phase_meas[deg] |Z_fit| Phase_fit[deg] ");
    fprintf(file, "Mag_Error[%%] Phase_Error[deg]\n");

    for (int i = 0; i < ctx->data->n_points; i++)
    {
        Complex z_fitted = ctx->model_func(ctx->data->frequencies[i], ctx->params);
        Complex z_measured = ctx->data->impedances[i];

        double mag_meas = sqrt(z_measured.real * z_measured.real +
                               z_measured.imag * z_measured.imag);
        double phase_meas = atan2(z_measured.imag, z_measured.real) * 180.0 / PI;

        double mag_fit = sqrt(z_fitted.real * z_fitted.real +
                              z_fitted.imag * z_fitted.imag);
        double phase_fit = atan2(z_fitted.imag, z_fitted.real) * 180.0 / PI;

        double mag_error_pct = 100.0 * (mag_fit - mag_meas) / mag_meas;
        double phase_error = phase_fit - phase_meas;
        while (phase_error > 180.0)
            phase_error -= 360.0;
        while (phase_error < -180.0)
            phase_error += 360.0;

        fprintf(file, "%.6e %.10e %.10e %.10e %.10e %.6e %.3f %.6e %.3f %.3f %.3f\n",
                ctx->data->frequencies[i],
                z_measured.real, z_measured.imag,
                z_fitted.real, z_fitted.imag,
                mag_meas, phase_meas, mag_fit, phase_fit,
                mag_error_pct, phase_error);
    }

    fclose(file);
    printf("Results exported to: %s\n", filename);
}

/* --- Print parameter sensitivity analysis --- */
void print_parameter_sensitivity(FittingContext *ctx, const double *par)
{
    if (!ctx || !par)
        return;

    printf("\nParameter Sensitivity Analysis:\n");
    print_separator('=', 50);

    const double delta = 0.01; /* 1% perturbation */
    double base_chi_sq = calculate_chi_squared(ctx, par);

    printf("%-10s %12s %12s %12s\n", "Parameter", "Value", "Chi2", "Sensitivity");
    print_separator('-', 50);

    for (int i = 0; i < ctx->params->n_params; i++)
    {
        double original_value = par[i];
        double *par_copy = malloc(ctx->params->n_params * sizeof(double));
        memcpy(par_copy, par, ctx->params->n_params * sizeof(double));

        /* Perturb parameter */
        par_copy[i] = original_value * (1.0 + delta);
        double perturbed_chi_sq = calculate_chi_squared(ctx, par_copy);

        double sensitivity = (perturbed_chi_sq - base_chi_sq) / (delta * original_value);

        printf("%-10s %12.6e %12.6e %12.6e\n",
               ctx->params->names[i], original_value, perturbed_chi_sq, sensitivity);

        free(par_copy);
    }
    printf("\n");
}