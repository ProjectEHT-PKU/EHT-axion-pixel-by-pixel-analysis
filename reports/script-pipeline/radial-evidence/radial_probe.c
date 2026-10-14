#include "Axion_Radial.h"
int main(void) {
    double spins[] = {0.9375, 0.99};
    double mus[] = {0.15, 0.3, 0.31, 0.4, 0.44};
    double radii[] = {1., 1.1, 1.2, 1.4, 2., 5., 10., 20., 100., 240.9};
    for (int i=0;i<2;i++) for(int j=0;j<5;j++) {
        double a=spins[i],mu=mus[j];
        if (a==0.9375 && mu>.316) continue;
        double complex w=get_omega(mu,a);
        printf("OMEGA %.4f %.4f %.17g %.17g\n",a,mu,creal(w),cimag(w));
        for(int k=0;k<10;k++) {
            double X[4]={0,radii[k],0,0};
            printf("R %.4f %.4f %.4f %.17g\n",a,mu,radii[k],Function_R(X,a,mu));
        }
    }
    return 0;
}
