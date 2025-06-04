#include "../src/eis_framework.h"

/* --- Simple Randles circuit model: Rs - (Rct || Cdl) --- */

typedef enum
{
    PARAM_RS = 0, /* Solution resistance */
    PARAM_RCT,    /* Charge transfer resistance */
    PARAM_CDL,    /* Double layer capacitance */
    N_RANDLES_PARAMS
} RandlesParamIndex;

/* --- Randles circuit model function --- */
Complex randles_model(double frequency, const Parameters *params)
{
    if (!params || params->n_params < N_RANDLES_PARAMS)
    {
        printf("Error: Invalid parameters for Randles model\n");
        return (Complex){0.0, 0.0};
    }

    double w = 2 * PI * frequency;

    double Rs = params->values[PARAM_RS];
    double Rct = params->values[PARAM_RCT];
    double Cdl = params->values[PARAM_CDL];

    Complex Z_Rs = resistor(Rs);
    Complex Z_Rct = resistor(Rct);
    Complex Z_Cdl = capacitor(Cdl, w);

    /* Rs in series with (Rct || Cdl) */
    Complex Z_rc = z_par(Z_Rct, Z_Cdl);
    return z_ser(Z_Rs, Z_rc);
}

/* --- Setup function for Randles model --- */
Parameters *setup_randles_parameters()
{
    Parameters *params = create_parameters(N_RANDLES_PARAMS);
    if (!params)
        return NULL;

    /* Set parameter properties with initial guesses of 1 */
    set_parameter(params, PARAM_RS, "Rs", 1.0, 1e-8, 1e8, LM_BOUND_NONE, 1000.0);
    set_parameter(params, PARAM_RCT, "Rct", 1.0, 1e-8, 1e8, LM_BOUND_NONE, 100.0);
    set_parameter(params, PARAM_CDL, "Cdl", 1.0, 1e-8, 1e8, LM_BOUND_NONE, 1e-6);

    return params;
}

/* --- Generate ideal Randles circuit data --- */
EISData *generate_randles_simulation_data()
{
    /* True parameter values */
    const double Rs_true = 3010.0; /* 3010 ohm */
    const double Rct_true = 100.0; /* 100 ohm */
    const double Cdl_true = 1e-6;  /* 1 uF */

    /* Frequency range: 0.01 Hz to 100 kHz (logarithmic) */
    const int n_points = 70;
    const double freq_min = 0.01;     /* 0.01 Hz */
    const double freq_max = 100000.0; /* 100 kHz */

    EISData *data = create_eis_data(n_points);
    if (!data)
        return NULL;

    /* Create true parameters for simulation */
    Parameters *true_params = create_parameters(N_RANDLES_PARAMS);
    if (!true_params)
    {
        free_eis_data(data);
        return NULL;
    }

    true_params->values[PARAM_RS] = Rs_true;
    true_params->values[PARAM_RCT] = Rct_true;
    true_params->values[PARAM_CDL] = Cdl_true;

    /* Generate logarithmic frequency range */
    double log_freq_min = log10(freq_min);
    double log_freq_max = log10(freq_max);

    printf("Generating Randles simulation data:\n");
    printf("True parameters: Rs = %.0f ohm, Rct = %.0f ohm, Cdl = %.0f uF\n",
           Rs_true, Rct_true, Cdl_true * 1e6);
    printf("Frequency range: %.2f Hz to %.0f Hz (%d points)\n",
           freq_min, freq_max, n_points);

    for (int i = 0; i < n_points; i++)
    {
        /* Logarithmic frequency spacing */
        double log_freq = log_freq_min + (log_freq_max - log_freq_min) * i / (n_points - 1);
        data->frequencies[i] = pow(10.0, log_freq);

        /* Calculate ideal impedance */
        data->impedances[i] = randles_model(data->frequencies[i], true_params);

        /* Add small amount of noise (0.1% of magnitude) */
        double noise_level = 0.001; /* 0.1% noise */
        double mag = sqrt(data->impedances[i].real * data->impedances[i].real +
                          data->impedances[i].imag * data->impedances[i].imag);

        /* Generate random noise */
        double noise_real = ((double)rand() / RAND_MAX - 0.5) * 2 * noise_level * mag;
        double noise_imag = ((double)rand() / RAND_MAX - 0.5) * 2 * noise_level * mag;

        data->impedances[i].real += noise_real;
        data->impedances[i].imag += noise_imag;

        /* Print first few points for verification */
        if (i < 5 || i >= n_points - 3)
        {
            printf("  f=%.3e Hz: Z = %.3e %+.3ej ohm\n",
                   data->frequencies[i],
                   data->impedances[i].real,
                   data->impedances[i].imag);
        }
        else if (i == 5)
        {
            printf("  ...\n");
        }
    }

    /* Set data source */
    data->data_source = malloc(50);
    if (data->data_source)
    {
        strcpy(data->data_source, "Simulated_Randles_Rs3010_Rct100_Cdl1uF");
    }

    free_parameters(true_params);
    return data;
}

/* --- Export simulation data to file --- */
void export_randles_simulation_data(const char *filename)
{
    EISData *data = generate_randles_simulation_data();
    if (!data)
    {
        printf("Failed to generate simulation data\n");
        return;
    }

    FILE *file = fopen(filename, "w");
    if (!file)
    {
        printf("Cannot create file: %s\n", filename);
        free_eis_data(data);
        return;
    }

    fprintf(file, "# Simulated Randles Circuit Data\n");
    fprintf(file, "# Circuit: Rs - (Rct || Cdl)\n");
    fprintf(file, "# True parameters: Rs=3010ohm, Rct=100ohm, Cdl=1uF\n");
    fprintf(file, "# Noise level: 0.1%% of magnitude\n");
    fprintf(file, "# Frequency[Hz] Z_real[ohm] Z_imag[ohm]\n");

    for (int i = 0; i < data->n_points; i++)
    {
        fprintf(file, "%.6e %.10e %.10e\n",
                data->frequencies[i],
                data->impedances[i].real,
                data->impedances[i].imag);
    }

    fclose(file);
    printf("Simulation data exported to: %s\n", filename);
    free_eis_data(data);
}

/* --- Fitting function for Randles model with known true values --- */
int fit_randles_to_simulation(const char *output_file)
{
    printf("Randles Circuit Fitting Test\n");
    printf("============================\n\n");

    /* Generate simulation data */
    EISData *data = generate_randles_simulation_data();
    if (!data)
    {
        printf("Failed to generate simulation data\n");
        return -1;
    }

    printf("\n");

    /* Setup parameters with initial guess */
    Parameters *params = setup_randles_parameters();
    if (!params)
    {
        free_eis_data(data);
        return -1;
    }

    /* Create fitting context */
    FittingContext *ctx = create_fitting_context(data, randles_model, params);
    if (!ctx)
    {
        free_parameters(params);
        free_eis_data(data);
        return -1;
    }

    /* Setup fitting control */
    lm_control_struct control = lm_control_double;
    control.stepbound = 0.1; /* Larger step bound for big parameter changes */
    control.patience = 2000; /* More iterations */
    control.scale_diag = 0;  /* Use parameter scaling */
    control.verbosity = 3;   /* Detailed output */
    control.ftol = 1e-12;
    control.xtol = 1e-12;
    control.gtol = 1e-12;

    /* Use modulus auto weighting */
    set_weighting_method(ctx, WEIGHT_UNITY);

    printf("Initial parameter guesses:\n");
    for (int i = 0; i < params->n_params; i++)
    {
        printf("  %-6s = %12.6g\n", params->names[i], params->values[i]);
    }
    printf("\n");

    printf("True parameter values:\n");
    printf("  Rs     = %12.0f ohm\n", 3010.0);
    printf("  Rct    = %12.0f ohm\n", 100.0);
    printf("  Cdl    = %12.0f uF\n", 1.0);
    printf("\n");

    /* Perform fitting */
    int result = fit_eis_model(ctx, &control);

    if (result >= 0)
    {
        printf("\nFitting Results vs True Values:\n");
        printf("===============================\n");

        double true_values[3] = {3010.0, 100.0, 1e-6};
        const char *units[3] = {"ohm", "ohm", "uF"};
        double unit_factors[3] = {1.0, 1.0, 1e6}; /* Convert Cdl to uF for display */

        printf("%-10s %12s %12s %12s %12s\n",
               "Parameter", "True", "Fitted", "Error%", "Unit");
        printf("------------------------------------------------------------\n");

        for (int i = 0; i < params->n_params; i++)
        {
            double fitted_display = params->values[i] * unit_factors[i];
            double true_display = true_values[i] * unit_factors[i];
            double error_pct = 100.0 * (params->values[i] - true_values[i]) / true_values[i];

            printf("%-10s %12.3f %12.3f %12.4f %12s\n",
                   params->names[i], true_display, fitted_display, error_pct, units[i]);
        }

        /* Check if fit is successful (within 1% error) */
        int fit_success = 1;
        for (int i = 0; i < params->n_params; i++)
        {
            double error_pct = fabs(100.0 * (params->values[i] - true_values[i]) / true_values[i]);
            if (error_pct > 1.0)
            {
                fit_success = 0;
                break;
            }
        }

        printf("\nFit quality assessment:\n");
        if (fit_success)
        {
            printf("EXCELLENT: All parameters within 1%% of true values\n");
        }
        else
        {
            printf("WARNING: Some parameters have >1%% error\n");
        }

        /* Show detailed analysis */
        print_frequency_comparison(ctx, params->values);

        /* Export results if requested */
        if (output_file)
        {
            export_fit_results(ctx, params->values, output_file);
            printf("Detailed results exported to: %s\n", output_file);
        }
    }

    /* Cleanup */
    free_fitting_context(ctx);
    free_parameters(params);
    free_eis_data(data);

    return result;
}

/* --- Test different initial guesses --- */
void test_initial_guess_sensitivity()
{
    printf("Testing Initial Guess Sensitivity\n");
    printf("=================================\n\n");

    /* Generate data once */
    EISData *data = generate_randles_simulation_data();
    if (!data)
        return;

    /* Test different initial guesses */
    double initial_guesses[][3] = {
        {1.0, 1.0, 1.0},         /* Original test case */
        {1000.0, 50.0, 1e-6},    /* Closer to true values */
        {10000.0, 1000.0, 1e-9}, /* Overestimate */
        {100.0, 10.0, 1e-3},     /* Underestimate */
        {0.1, 0.1, 1e-12}        /* Very low initial guess */
    };

    const char *guess_names[] = {
        "All ones", "Reasonable guess", "Overestimate", "Underestimate", "Very low"};

    int n_tests = sizeof(initial_guesses) / sizeof(initial_guesses[0]);

    for (int test = 0; test < n_tests; test++)
    {
        printf("Test %d: %s\n", test + 1, guess_names[test]);
        printf("Initial: Rs=%.1e, Rct=%.1e, Cdl=%.1e\n",
               initial_guesses[test][0], initial_guesses[test][1], initial_guesses[test][2]);

        Parameters *params = setup_randles_parameters();
        params->values[PARAM_RS] = initial_guesses[test][0];
        params->values[PARAM_RCT] = initial_guesses[test][1];
        params->values[PARAM_CDL] = initial_guesses[test][2];

        FittingContext *ctx = create_fitting_context(data, randles_model, params);
        set_weighting_method(ctx, WEIGHT_UNITY);

        lm_control_struct control = lm_control_double;
        control.stepbound = 0.1;
        control.patience = 1000;
        control.scale_diag = 0;
        control.verbosity = 1; /* Less verbose for multiple tests */
        control.ftol = 1e-12;
        control.xtol = 1e-12;
        control.gtol = 1e-12;

        int result = fit_eis_model(ctx, &control);

        if (result >= 0)
        {
            double true_values[3] = {3010.0, 100.0, 1e-6};
            double max_error = 0.0;

            for (int i = 0; i < 3; i++)
            {
                double error_pct = fabs(100.0 * (params->values[i] - true_values[i]) / true_values[i]);
                if (error_pct > max_error)
                    max_error = error_pct;
            }

            printf("Result: Rs=%.1f, Rct=%.1f, Cdl=%.1f uF\n",
                   params->values[0], params->values[1], params->values[2] * 1e6);
            printf("Max error: %.3f%%\n", max_error);

            if (max_error < 1.0)
            {
                printf("Status: SUCCESS\n");
            }
            else
            {
                printf("Status: FAILED\n");
            }
        }
        else
        {
            printf("Status: OPTIMIZATION FAILED\n");
        }

        printf("\n");

        free_fitting_context(ctx);
        free_parameters(params);
    }

    free_eis_data(data);
}