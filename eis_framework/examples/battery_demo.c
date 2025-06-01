#include "../src/eis_framework.h"

/* External functions from battery_model.c */
extern Complex battery_model(double frequency, const Parameters *params);
extern Parameters *setup_battery_parameters();
extern EISData *load_battery_data();
extern void export_battery_simulation_data(const char *filename);
extern int fit_battery_model(const char *output_file);
extern void compare_battery_with_randles();

/* --- Circuit diagram for original model --- */
void print_battery_circuit()
{
    printf("EIS Model Circuit\n");
    printf("==================================\n\n");
    printf("Circuit: Rstray||L1 - Rs - Cf||Rf - Cdl||(Rct+Q1)\n\n");
    printf("Component descriptions:\n");
    printf("  Rstray = Stray resistance (parasitic)\n");
    printf("  L1     = Inductance (leads, connections)\n");
    printf("  Rs     = Series resistance (electrolyte)\n");
    printf("  Cf     = Film capacitance\n");
    printf("  Rf     = Film resistance\n");
    printf("  Cdl    = Double layer capacitance\n");
    printf("  Rct    = Charge transfer resistance\n");
    printf("  Q1     = Constant Phase Element (CPE)\n\n");
}

/* --- Interactive menu --- */
void show_battery_menu()
{
    printf("EIS Model Menu:\n");
    printf("=======================\n\n");
    printf("1. Export original eis.c data to file\n");
    printf("2. Fit original model (reproduce eis.c results)\n");
    printf("3. Show circuit diagram\n");
    printf("4. Compare with simple Randles model\n");
    printf("5. Run complete analysis\n");
    printf("0. Exit\n\n");
    printf("Choose option: ");
}

int main(int argc, char *argv[])
{
    printf("Original EIS Model Demonstration (eis.c reproduction)\n");
    printf("====================================================\n\n");

    if (argc > 1)
    {
        /* Command line mode */
        if (strcmp(argv[1], "export") == 0)
        {
            const char *filename = (argc > 2) ? argv[2] : "original_eis_data.txt";
            export_battery_simulation_data(filename);
        }
        else if (strcmp(argv[1], "fit") == 0)
        {
            const char *output = (argc > 2) ? argv[2] : "original_fit_results.txt";
            fit_battery_model(output);
        }
        else if (strcmp(argv[1], "compare") == 0)
        {
            compare_battery_with_randles();
        }
        else if (strcmp(argv[1], "all") == 0)
        {
            print_battery_circuit();
            export_battery_simulation_data("original_eis_data.txt");
            fit_battery_model("original_fit_results.txt");
            compare_battery_with_randles();
        }
        else
        {
            printf("Usage: %s [export|fit|compare|all] [output_file]\n", argv[0]);
            return 1;
        }
    }
    else
    {
        /* Interactive mode */
        int choice;
        do
        {
            show_battery_menu();
            if (scanf("%d", &choice) != 1)
            {
                printf("Invalid input!\n");
                while (getchar() != '\n')
                    ;
                continue;
            }

            printf("\n");
            switch (choice)
            {
            case 1:
                export_battery_simulation_data("original_eis_data.txt");
                break;
            case 2:
                fit_battery_model("original_fit_results.txt");
                break;
            case 3:
                print_battery_circuit();
                break;
            case 4:
                compare_battery_with_randles();
                break;
            case 5:
                print_battery_circuit();
                export_battery_simulation_data("original_eis_data.txt");
                fit_battery_model("original_fit_results.txt");
                compare_battery_with_randles();
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
                    ;
                getchar();
                printf("\n");
            }

        } while (choice != 0);
    }

    return 0;
}