#include "../src/eis_framework.h"

void run_simple_demo()
{
    printf("EIS Framework Simple Demo\n");
    printf("========================\n\n");

    printf("Available demos:\n");
    printf("  eis_randles_demo.exe  - Simple Randles circuit\n");
    printf("  eis_battery_demo.exe  - Original eis.c model\n\n");

    printf("Example usage:\n");
    printf("  eis_randles_demo.exe fit\n");
    printf("  eis_battery_demo.exe fit\n");
    printf("  eis_battery_demo.exe compare\n\n");
}

int main(int argc, char *argv[])
{
    printf("EIS Framework Main Example\n");
    printf("=========================\n\n");

    if (argc > 1 && strcmp(argv[1], "help") == 0)
    {
        printf("This is a placeholder main example.\n");
        printf("Please use the specific demos instead:\n\n");
    }

    run_simple_demo();

    return 0;
}