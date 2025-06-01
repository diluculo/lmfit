#include "../src/eis_framework.h"

/* External functions from randles_model.c */
extern Complex randles_model(double frequency, const Parameters *params);
extern Parameters *setup_randles_parameters();
extern EISData *generate_randles_simulation_data();
extern void export_randles_simulation_data(const char *filename);
extern int fit_randles_to_simulation(const char *output_file);
extern void test_initial_guess_sensitivity();

/* --- Circuit diagram for Randles model --- */
void print_randles_circuit()
{
    printf("Simple Randles Circuit\n");
    printf("Rs  = Solution resistance\n");
    printf("Rct = Charge transfer resistance\n");
    printf("Cdl = Double layer capacitance\n\n");
}

/* --- Nyquist plot prediction --- */
void predict_nyquist_behavior()
{
    printf("Expected Nyquist Plot Behavior:\n");
    printf("==============================\n\n");

    printf("For Rs=3010Ω, Rct=100Ω, Cdl=1µF:\n\n");

    printf("High frequency limit (f→∞):\n");
    printf("  - Capacitor acts as short circuit\n");
    printf("  - Z → Rs = 3010 Ω\n");
    printf("  - Real axis intercept at left\n\n");

    printf("Low frequency limit (f→0):\n");
    printf("  - Capacitor acts as open circuit\n");
    printf("  - Z → Rs + Rct = 3110 Ω\n");
    printf("  - Real axis intercept at right\n\n");

    printf("Characteristic frequency:\n");
    printf("  fc = 1/(2π·Rct·Cdl) = 1/(2π·100·1e-6) ≈ 1592 Hz\n");
    printf("  - Maximum imaginary component at this frequency\n");
    printf("  - Semicircle center: (3060, 0)\n");
    printf("  - Semicircle radius: 50 Ω\n\n");
}

/* --- Interactive menu --- */
void show_menu()
{
    printf("Randles Model Demonstration Menu:\n");
    printf("================================\n\n");
    printf("1. Generate and export simulation data\n");
    printf("2. Fit simulation data (starting from Rs=1, Rct=1, Cdl=1)\n");
    printf("3. Test initial guess sensitivity\n");
    printf("4. Show circuit diagram and theory\n");
    printf("5. Predict Nyquist plot behavior\n");
    printf("6. Run complete analysis\n");
    printf("0. Exit\n\n");
    printf("Choose option: ");
}

int main(int argc, char *argv[])
{
    printf("Randles Circuit EIS Fitting Demonstration\n");
    printf("=========================================\n\n");

    /* Seed random number generator for consistent noise */
    srand(12345);

    if (argc > 1)
    {
        /* Command line mode */
        if (strcmp(argv[1], "generate") == 0)
        {
            const char *filename = (argc > 2) ? argv[2] : "randles_simulation.txt";
            export_randles_simulation_data(filename);
        }
        else if (strcmp(argv[1], "fit") == 0)
        {
            const char *output = (argc > 2) ? argv[2] : "randles_fit_results.txt";
            fit_randles_to_simulation(output);
        }
        else if (strcmp(argv[1], "test") == 0)
        {
            test_initial_guess_sensitivity();
        }
        else if (strcmp(argv[1], "all") == 0)
        {
            print_randles_circuit();
            predict_nyquist_behavior();
            export_randles_simulation_data("randles_simulation.txt");
            fit_randles_to_simulation("randles_fit_results.txt");
            test_initial_guess_sensitivity();
        }
        else
        {
            printf("Usage: %s [generate|fit|test|all] [output_file]\n", argv[0]);
            return 1;
        }
    }
    else
    {
        /* Interactive mode */
        int choice;
        do
        {
            show_menu();
            if (scanf("%d", &choice) != 1)
            {
                printf("Invalid input!\n");
                while (getchar() != '\n')
                    ; /* Clear input buffer */
                continue;
            }

            printf("\n");
            switch (choice)
            {
            case 1:
                export_randles_simulation_data("randles_simulation.txt");
                break;
            case 2:
                fit_randles_to_simulation("randles_fit_results.txt");
                break;
            case 3:
                test_initial_guess_sensitivity();
                break;
            case 4:
                print_randles_circuit();
                break;
            case 5:
                predict_nyquist_behavior();
                break;
            case 6:
                print_randles_circuit();
                predict_nyquist_behavior();
                export_randles_simulation_data("randles_simulation.txt");
                fit_randles_to_simulation("randles_fit_results.txt");
                test_initial_guess_sensitivity();
                break;
            case 0:
                printf("Goodbye!\n");
                break;
            default:
                printf("Invalid choice!\n");
            }

            if (choice != 0)
            {
                printf("\nPress Enter to continue...");
                while (getchar() != '\n')
                    ;      /* Clear buffer */
                getchar(); /* Wait for Enter */
                printf("\n");
            }

        } while (choice != 0);
    }

    return 0;
}