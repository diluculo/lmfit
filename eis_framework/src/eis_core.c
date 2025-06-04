#include "eis_framework.h"

/* --- Complex arithmetic functions --- */
Complex c_add(Complex a, Complex b)
{
    return (Complex){a.real + b.real, a.imag + b.imag};
}

Complex c_sub(Complex a, Complex b)
{
    return (Complex){a.real - b.real, a.imag - b.imag};
}

Complex c_mul(Complex a, Complex b)
{
    return (Complex){
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real};
}

Complex c_mul_scalar(Complex a, double b)
{
    return (Complex){a.real * b, a.imag * b};
}

Complex c_div(Complex a, Complex b)
{
    double denom = b.real * b.real + b.imag * b.imag;
    if (denom < 1e-20)
    {
        printf("Warning: Division by very small complex number\n");
        return (Complex){INFINITY, INFINITY};
    }
    return (Complex){
        (a.real * b.real + a.imag * b.imag) / denom,
        (a.imag * b.real - a.real * b.imag) / denom};
}

Complex c_inv(Complex a)
{
    double denom = a.real * a.real + a.imag * a.imag;
    if (denom < 1e-20)
    {
        printf("Warning: Inversion of very small complex number\n");
        return (Complex){INFINITY, INFINITY};
    }
    return (Complex){a.real / denom, -a.imag / denom};
}

/* --- Impedance composition functions --- */
Complex z_ser(Complex z1, Complex z2)
{
    return c_add(z1, z2);
}

Complex z_ser_array(int count, const Complex *z)
{
    Complex sum = {0.0, 0.0};
    for (int i = 0; i < count; ++i)
    {
        sum = z_ser(sum, z[i]);
    }
    return sum;
}

Complex z_par(Complex z1, Complex z2)
{
    Complex numerator = c_mul(z1, z2);
    Complex denominator = c_add(z1, z2);
    if (fabs(denominator.real) < 1e-20 && fabs(denominator.imag) < 1e-20)
    {
        printf("Warning: Parallel impedance with very small denominator\n");
        return (Complex){INFINITY, INFINITY};
    }
    return c_div(numerator, denominator);
}

/* --- Component model functions --- */
Complex resistor(double R)
{
    return (Complex){R, 0.0};
}

Complex inductor(double L, double omega)
{
    return (Complex){0.0, omega * L};
}

Complex capacitor(double C, double omega)
{
    return (Complex){0.0, -1.0 / (omega * C)};
}

Complex constantPhaseElement(double Q, double alpha, double omega)
{
    if (fabs(Q) < 1e-20 || fabs(omega) < 1e-20)
    {
        return (Complex){INFINITY, INFINITY};
    }

    double magnitude = pow(omega, alpha);
    double phase = alpha * PI / 2.0;

    Complex jw_alpha = {magnitude * cos(phase), magnitude * sin(phase)};
    Complex qz = c_mul_scalar(jw_alpha, Q);
    return c_inv(qz);
}

Complex warburg(double sigma, double omega)
{
    if (fabs(omega) < 1e-20)
    {
        return (Complex){INFINITY, INFINITY};
    }

    double sqrt_omega = sqrt(omega);
    double factor = sigma / sqrt_omega;
    return (Complex){factor, -factor};
}

/* --- Data management functions --- */
EISData *create_eis_data(int n_points)
{
    if (n_points <= 0 || n_points > MAX_FREQUENCIES)
    {
        printf("Invalid number of data points: %d\n", n_points);
        return NULL;
    }

    EISData *data = malloc(sizeof(EISData));
    if (!data)
    {
        printf("Failed to allocate EISData structure\n");
        return NULL;
    }

    data->n_points = n_points;
    data->frequencies = malloc(n_points * sizeof(double));
    data->impedances = malloc(n_points * sizeof(Complex));
    data->data_source = NULL;

    if (!data->frequencies || !data->impedances)
    {
        printf("Failed to allocate data arrays\n");
        free_eis_data(data);
        return NULL;
    }

    return data;
}

void free_eis_data(EISData *data)
{
    if (data)
    {
        free(data->frequencies);
        free(data->impedances);
        free(data->data_source);
        free(data);
    }
}

int load_data_from_arrays(double *freq, Complex *imp, int n_points, EISData *data)
{
    if (!data || !freq || !imp || n_points != data->n_points)
    {
        return -1;
    }

    memcpy(data->frequencies, freq, n_points * sizeof(double));
    memcpy(data->impedances, imp, n_points * sizeof(Complex));

    return 0;
}

int load_data_from_file(const char *filename, EISData *data)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        printf("Cannot open file: %s\n", filename);
        return -1;
    }

    // Count lines first
    int n_lines = 0;
    char line[1024];
    while (fgets(line, sizeof(line), file))
    {
        if (line[0] != '#' && strlen(line) > 1)
        {
            n_lines++;
        }
    }

    if (n_lines == 0)
    {
        printf("No data found in file: %s\n", filename);
        fclose(file);
        return -1;
    }

    // Reallocate if necessary
    if (data->n_points != n_lines)
    {
        free(data->frequencies);
        free(data->impedances);
        data->n_points = n_lines;
        data->frequencies = malloc(n_lines * sizeof(double));
        data->impedances = malloc(n_lines * sizeof(Complex));
        if (!data->frequencies || !data->impedances)
        {
            printf("Failed to allocate memory for data\n");
            fclose(file);
            return -1;
        }
    }

    // Read data
    rewind(file);
    int i = 0;
    while (fgets(line, sizeof(line), file) && i < n_lines)
    {
        if (line[0] == '#' || strlen(line) <= 1)
            continue;

        double freq, real, imag;
        if (sscanf(line, "%lf %lf %lf", &freq, &real, &imag) == 3)
        {
            data->frequencies[i] = freq;
            data->impedances[i].real = real;
            data->impedances[i].imag = imag;
            i++;
        }
    }

    fclose(file);

    // Set data source
    if (data->data_source)
        free(data->data_source);
    data->data_source = malloc(strlen(filename) + 1);
    if (data->data_source)
    {
        strcpy(data->data_source, filename);
    }

    return 0;
}

/* --- Parameter management functions --- */
Parameters *create_parameters(int n_params)
{
    if (n_params <= 0 || n_params > MAX_PARAMS)
    {
        printf("Invalid number of parameters: %d\n", n_params);
        return NULL;
    }

    Parameters *params = malloc(sizeof(Parameters));
    if (!params)
    {
        printf("Failed to allocate Parameters structure\n");
        return NULL;
    }

    params->n_params = n_params;
    params->values = malloc(n_params * sizeof(double));
    params->names = malloc(n_params * sizeof(char *));
    params->lower_bounds = malloc(n_params * sizeof(double));
    params->upper_bounds = malloc(n_params * sizeof(double));
    params->bound_types = malloc(n_params * sizeof(int));
    params->scales = malloc(n_params * sizeof(double));

    if (!params->values || !params->names || !params->lower_bounds ||
        !params->upper_bounds || !params->bound_types || !params->scales)
    {
        printf("Failed to allocate parameter arrays\n");
        free_parameters(params);
        return NULL;
    }

    // Initialize names to NULL
    for (int i = 0; i < n_params; i++)
    {
        params->names[i] = NULL;
    }

    return params;
}

void free_parameters(Parameters *params)
{
    if (params)
    {
        free(params->values);
        if (params->names)
        {
            for (int i = 0; i < params->n_params; i++)
            {
                free(params->names[i]);
            }
            free(params->names);
        }
        free(params->lower_bounds);
        free(params->upper_bounds);
        free(params->bound_types);
        free(params->scales);
        free(params);
    }
}

int set_parameter(Parameters *params, int index, const char *name, double value,
                  double lower, double upper, int bound_type, double scale)
{
    if (!params || index < 0 || index >= params->n_params)
    {
        return -1;
    }

    params->values[index] = value;
    params->lower_bounds[index] = lower;
    params->upper_bounds[index] = upper;
    params->bound_types[index] = bound_type;
    params->scales[index] = scale;

    // Set name
    if (params->names[index])
    {
        free(params->names[index]);
    }
    params->names[index] = malloc(strlen(name) + 1);
    if (params->names[index])
    {
        strcpy(params->names[index], name);
    }

    return 0;
}