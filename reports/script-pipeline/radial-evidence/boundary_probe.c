#include "Axion_Radial.h"
int main(void) {
 double complex w=get_omega(.1,.9375);
 printf("%.17g %.17g\n",creal(w),cimag(w));
 return 0;
}
