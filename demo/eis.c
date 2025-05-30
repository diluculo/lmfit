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

void evaluate_residuals(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
    (void)info;
    (void)data;

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
        Complex z = model_impedance(frequencies[i], &p);
        fvec[i] = z.real - impedances[i].real;
        fvec[i + N] = z.imag - impedances[i].imag;
    }
}

/* --- Main function --- */

int main()
{
    const int n_par = 9;

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

    lm_control_struct control = lm_control_double;
    control.stepbound = 0.01; // Initial step bound
    control.patience = 200;   // Maximum number of function evaluations
    control.scale_diag = 1;   // Rescale variables internally

    lm_status_struct status;
    control.verbosity = 2;

    lmmin(n_par, par, 2 * N, y, NULL, evaluate_residuals, &control, &status);

    printf("\nFitted parameters:\n");
    for (int i = 0; i < n_par; ++i)
        printf("  par[%d] = %.10g\n", i, par[i]);

    return 0;
}
