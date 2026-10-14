#include "Axion_Radial.h"
int main(void) {
    double complex omega = get_omega(0.1, 0.9375);
    printf("omega_real=%.17g omega_imag=%.17g\n", creal(omega), cimag(omega));
    print_cable(0.9375, 0.1, 0.1);
    return 0;
}
