#include "../src/eis_framework.h"

/* --- Original EIS model from eis.c --- */
/* Circuit: Rstray|L1-Rs-Cf|Rf-Cdl|(Rct-Q1) */

typedef enum
{
    PARAM_RSTRAY = 0, /* Stray resistance */
    PARAM_L1,         /* Inductance */
    PARAM_RS,         /* Series resistance */
    PARAM_CF,         /* Film capacitance */
    PARAM_RF,         /* Film resistance */
    PARAM_CDL,        /* Double layer capacitance */
    PARAM_RCT,        /* Charge transfer resistance */
    PARAM_QY1,        /* CPE magnitude */
    PARAM_QA1,        /* CPE exponent */
    N_BATTERY_PARAMS
} BatteryParamIndex;

/* --- Battery model impedance function (original eis.c circuit) --- */
Complex battery_model(double frequency, const Parameters *params)
{
    if (!params || params->n_params < N_BATTERY_PARAMS)
    {
        printf("Error: Invalid parameters for battery model\n");
        return (Complex){0.0, 0.0};
    }

    double w = 2 * PI * frequency;

    double Rstray = params->values[PARAM_RSTRAY];
    double L1 = params->values[PARAM_L1];
    double Rs = params->values[PARAM_RS];
    double Cf = params->values[PARAM_CF];
    double Rf = params->values[PARAM_RF];
    double Cdl = params->values[PARAM_CDL];
    double Rct = params->values[PARAM_RCT];
    double Qy1 = params->values[PARAM_QY1];
    double Qa1 = params->values[PARAM_QA1];

    /* Build circuit components */
    Complex Z_Rstray = resistor(Rstray);
    Complex Z_L1 = inductor(L1, w);
    Complex Z_Rs = resistor(Rs);
    Complex Z_Cf = capacitor(Cf, w);
    Complex Z_Rf = resistor(Rf);
    Complex Z_Cdl = capacitor(Cdl, w);
    Complex Z_Rct = resistor(Rct);
    Complex Z_Q1 = constantPhaseElement(Qy1, Qa1, w);

    /* Circuit topology: Rstray|L1-Rs-Cf|Rf-Cdl|(Rct-Q1) */

    /* Z1 = Rstray || L1 */
    Complex Z1 = z_par(Z_Rstray, Z_L1);

    /* Z2 = Cf || Rf */
    Complex Z2 = z_par(Z_Cf, Z_Rf);

    /* Z3 = Cdl || (Rct + Q1) */
    Complex Z_series_rct_q1 = z_ser(Z_Rct, Z_Q1);
    Complex Z3 = z_par(Z_Cdl, Z_series_rct_q1);

    /* Total impedance: Z1 + Rs + Z2 + Z3 */
    Complex elements[4] = {Z1, Z_Rs, Z2, Z3};
    return z_ser_array(4, elements);
}

/* --- Setup function for battery model parameters --- */
Parameters *setup_battery_parameters()
{
    Parameters *params = create_parameters(N_BATTERY_PARAMS);
    if (!params)
        return NULL;

    /* Set parameter properties with initial values from original code */
    set_parameter(params, PARAM_RSTRAY, "Rstray", 10e-3, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1e-3);
    set_parameter(params, PARAM_L1, "L1", 103.7e-9, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1e-6);
    set_parameter(params, PARAM_RS, "Rs", 408.6e-6, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1e-3);
    set_parameter(params, PARAM_CF, "Cf", 109.4, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1.0);
    set_parameter(params, PARAM_RF, "Rf", 64.7e-6, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1e-3);
    set_parameter(params, PARAM_CDL, "Cdl", 6.185, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1.0);
    set_parameter(params, PARAM_RCT, "Rct", 171.4e-6, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1e-3);
    set_parameter(params, PARAM_QY1, "Qy1", 2500, 0.0, INFINITY,
                  LM_BOUND_LOWER, 1000.0);
    set_parameter(params, PARAM_QA1, "Qa1", 0.7, 0.0, 1.0,
                  LM_BOUND_BOTH, 0.1);

    return params;
}

/* --- Load original measured data --- */
EISData *load_battery_data()
{
    /* Original frequencies from eis.c */
    static const double frequencies[29] = {
        4000, 3162.28, 2154.43, 1467.8, 1000, 681.292, 464.159, 316.228,
        215.443, 146.78, 100, 68.129, 46.416, 31.623, 21.544, 14.678,
        10, 6.813, 4.641, 3.162, 2.154, 1.468, 1, 0.681,
        0.464, 0.316, 0.215, 0.147, 0.1};

    /* Original impedances from eis.c */
    static const Complex impedances[29] = {
        {0.000663089074, 0.00252012886}, {0.000593130962, 0.002012695201}, {0.000514193238, 0.001402412198}, {0.000451633546, 0.000951764771}, {0.000425432787, 0.00062922476}, {0.000420608933, 0.000410356128}, {0.000427075089, 0.000249392344}, {0.000443342498, 0.000133066927}, {0.000466690484, 0.000052738781}, {0.000496198248, 0.000000873319}, {0.000523345566, -0.000028826491}, {0.000551174431, -0.000043343169}, {0.000580768587, -0.000054246608}, {0.000601011623, -0.000057778532}, {0.000616794061, -0.000058350616}, {0.000628636927, -0.000055934551}, {0.000641708479, -0.000056641414}, {0.000652100449, -0.000061967468}, {0.000662929832, -0.000066026}, {0.00067382404, -0.000075871667}, {0.000689054979, -0.0000893174}, {0.000702894694, -0.000108925309}, {0.000719473295, -0.000131383953}, {0.00073985794, -0.000171072749}, {0.000766850441, -0.00021927308}, {0.000801383826, -0.00028480345}, {0.000845546901, -0.000375287209}, {0.000914403545, -0.000494018774}, {0.001018981637, -0.000651889399}};

    const int n_points = 29;

    EISData *data = create_eis_data(n_points);
    if (!data)
        return NULL;

    /* Copy frequency and impedance data */
    for (int i = 0; i < n_points; i++)
    {
        data->frequencies[i] = frequencies[i];
        data->impedances[i] = impedances[i];
    }

    /* Set data source */
    data->data_source = malloc(30);
    if (data->data_source)
    {
        strcpy(data->data_source, "Original_eis.c_dataset");
    }

    return data;
}

/* --- Export original data to file --- */
void export_battery_simulation_data(const char *filename)
{
    EISData *data = load_battery_data();
    if (!data)
    {
        printf("Failed to load original data\n");
        return;
    }

    FILE *file = fopen(filename, "w");
    if (!file)
    {
        printf("Cannot create file: %s\n", filename);
        free_eis_data(data);
        return;
    }

    fprintf(file, "# Original EIS data from eis.c\n");
    fprintf(file, "# Circuit: Rstray||L1 - Rs - Cf||Rf - Cdl||(Rct+Q1)\n");
    fprintf(file, "# Complex electrochemical system\n");
    fprintf(file, "# 29 frequency points from 4000 Hz to 0.1 Hz\n");
    fprintf(file, "# Frequency[Hz] Z_real[Ω] Z_imag[Ω]\n");

    for (int i = 0; i < data->n_points; i++)
    {
        fprintf(file, "%.6e %.12e %.12e\n",
                data->frequencies[i],
                data->impedances[i].real,
                data->impedances[i].imag);
    }

    fclose(file);
    printf("Original EIS data exported to: %s\n", filename);
    free_eis_data(data);
}

/* --- Fitting function for battery model --- */
int fit_battery_model(const char *output_file)
{
    printf("Original EIS Model Fitting (eis.c reproduction)\n");
    printf("===============================================\n\n");

    /* Load original data */
    EISData *data = load_battery_data();
    if (!data)
    {
        printf("Failed to load original data\n");
        return -1;
    }

    printf("Loaded %d data points from original eis.c\n", data->n_points);
    printf("Frequency range: %.3f Hz to %.0f Hz\n",
           data->frequencies[data->n_points - 1], data->frequencies[0]);
    printf("Circuit: Rstray||L1 - Rs - Cf||Rf - Cdl||(Rct+Q1)\n\n");

    /* Setup parameters */
    Parameters *params = setup_battery_parameters();
    if (!params)
    {
        free_eis_data(data);
        return -1;
    }

    /* Create fitting context */
    FittingContext *ctx = create_fitting_context(data, battery_model, params);
    if (!ctx)
    {
        free_parameters(params);
        free_eis_data(data);
        return -1;
    }

    /* Setup fitting control - same as original */
    lm_control_struct control = lm_control_double;
    control.stepbound = 0.1;
    control.patience = 1000;
    control.scale_diag = 1;
    control.verbosity = 3;
    control.ftol = 1e-15;
    control.xtol = 1e-15;
    control.gtol = 1e-15;

    /* Use modulus auto weighting (same as original) */
    set_weighting_method(ctx, WEIGHT_MODULUS_AUTO);

    printf("Initial parameter values:\n");
    for (int i = 0; i < params->n_params; i++)
    {
        printf("  %-6s = %12.6g\n", params->names[i], params->values[i]);
    }
    printf("\n");

    /* Perform fitting */
    int result = fit_eis_model(ctx, &control);

    if (result >= 0)
    {
        /* Show detailed analysis */
        print_frequency_comparison(ctx, params->values);
        print_parameter_sensitivity(ctx, params->values);

        /* Export results if requested */
        if (output_file)
        {
            export_fit_results(ctx, params->values, output_file);
        }

        /* Calculate total squared error for comparison with original */
        double total_error = 0.0;
        for (int i = 0; i < data->n_points; i++)
        {
            /* Update parameters in case they changed during fitting */
            for (int j = 0; j < params->n_params; j++)
            {
                params->values[j] = params->values[j];
            }

            Complex z_model = battery_model(data->frequencies[i], params);
            Complex z_meas = data->impedances[i];

            double error_real = z_model.real - z_meas.real;
            double error_imag = z_model.imag - z_meas.imag;
            total_error += error_real * error_real + error_imag * error_imag;
        }

        printf("\nComparison with original eis.c:\n");
        printf("==============================\n");
        printf("Total squared error: %.6e\n", total_error);
        printf("(Should match original eis.c output)\n\n");

        printf("Expected final parameter values from original eis.c:\n");
        printf("(Compare with fitted values above)\n");
        printf("  Rstray ≈ 10e-3\n");
        printf("  L1     ≈ 103.7e-9\n");
        printf("  Rs     ≈ 408.6e-6\n");
        printf("  Cf     ≈ 109.4\n");
        printf("  Rf     ≈ 64.7e-6\n");
        printf("  Cdl    ≈ 6.185\n");
        printf("  Rct    ≈ 171.4e-6\n");
        printf("  Qy1    ≈ 2500\n");
        printf("  Qa1    ≈ 0.7\n");
    }

    /* Cleanup */
    free_fitting_context(ctx);
    free_parameters(params);
    free_eis_data(data);

    return result;
}

/* --- Function to compare with simpler models --- */
void compare_battery_with_randles()
{
    printf("Comparing Original Complex Model with Simple Randles\n");
    printf("====================================================\n\n");

    EISData *data = load_battery_data();
    if (!data)
        return;

    /* Fit original complex model */
    Parameters *battery_params = setup_battery_parameters();
    FittingContext *battery_ctx = create_fitting_context(data, battery_model, battery_params);

    lm_control_struct control = lm_control_double;
    control.stepbound = 0.1;
    control.patience = 1000;
    control.scale_diag = 1;
    control.verbosity = 1;
    control.ftol = 1e-15;
    control.xtol = 1e-15;
    control.gtol = 1e-15;

    set_weighting_method(battery_ctx, WEIGHT_MODULUS_AUTO);

    printf("Fitting original complex model (9 parameters)...\n");
    fit_eis_model(battery_ctx, &control);
    double battery_chi_sq = calculate_chi_squared(battery_ctx, battery_params->values);

    /* Fit simple Randles model for comparison */
    extern Parameters *setup_randles_parameters(); /* From randles_model.c */
    extern Complex randles_model(double frequency, const Parameters *params);

    Parameters *randles_params = setup_randles_parameters();
    FittingContext *randles_ctx = create_fitting_context(data, randles_model, randles_params);
    set_weighting_method(randles_ctx, WEIGHT_MODULUS_AUTO);

    printf("\nFitting simple Randles model (3 parameters)...\n");
    fit_eis_model(randles_ctx, &control);
    double randles_chi_sq = calculate_chi_squared(randles_ctx, randles_params->values);

    /* Compare results */
    printf("\nModel Comparison Results:\n");
    printf("========================\n");
    printf("Original complex model (9 params): χ² = %.6e\n", battery_chi_sq);
    printf("Simple Randles model   (3 params): χ² = %.6e\n", randles_chi_sq);
    printf("χ² ratio (Randles/Original): %.2f\n", randles_chi_sq / battery_chi_sq);

    if (randles_chi_sq / battery_chi_sq < 2.0)
    {
        printf("→ Simple Randles model provides comparable fit!\n");
    }
    else
    {
        printf("→ Original complex model provides significantly better fit.\n");
        printf("→ The additional parameters are justified by the data.\n");
    }

    /* Cleanup */
    free_fitting_context(battery_ctx);
    free_fitting_context(randles_ctx);
    free_parameters(battery_params);
    free_parameters(randles_params);
    free_eis_data(data);
}