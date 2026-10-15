#include <stdio.h>
#include <math.h>
#include <complex.h>
int check_stiffness(double jI, double jQ, double jU, double jV, double rQ,
                    double rU, double rV, double aI, double aQ, double aU,
                    double aV, double dl_current) {
    // int STIFF = check_stiffness...
    double a2 = rQ * rQ + rV * rV - aQ * aQ - aV * aV;
    double a0 = -2. * aV * aQ * rV * rQ - aQ * aQ * rQ * rQ - aV * aV * rV * rV;

    complex double zplus = (-a2 + sqrt(a2 * a2 - 4. * a0)) / 2.;
    complex double zminus = (-a2 - sqrt(a2 * a2 - 4. * a0)) / 2.;

    complex double l1 = aI + sqrt(zplus);
    complex double l2 = aI - sqrt(zplus);
    complex double l3 = aI + sqrt(zminus);
    complex double l4 = aI - sqrt(zminus);

    complex double tau1 = dl_current * l1;
    complex double tau2 = dl_current * l2;
    complex double tau3 = dl_current * l3;
    complex double tau4 = dl_current * l4;

    complex double mag1 = 1. + tau1 + 0.5 * tau1 * tau1 +
                          1. / 6. * tau1 * tau1 * tau1 +
                          1. / 24. * tau1 * tau1 * tau1 * tau1;
    complex double mag2 = 1. + tau2 + 0.5 * tau2 * tau2 +
                          1. / 6. * tau2 * tau2 * tau2 +
                          1. / 24. * tau2 * tau2 * tau2 * tau2;
    complex double mag3 = 1. + tau3 + 0.5 * tau3 * tau3 +
                          1. / 6. * tau3 * tau3 * tau3 +
                          1. / 24. * tau3 * tau3 * tau3 * tau3;
    complex double mag4 = 1. + tau4 + 0.5 * tau4 * tau4 +
                          1. / 6. * tau4 * tau4 * tau4 +
                          1. / 24. * tau4 * tau4 * tau4 * tau4;

    double res1 = sqrt(mag1 * conj(mag1));
    double res2 = sqrt(mag2 * conj(mag2));
    double res3 = sqrt(mag3 * conj(mag3));
    double res4 = sqrt(mag4 * conj(mag4));

    int STIFF = 0;

    double STIFFTHRESH = 0.99;

    if (res1 > STIFFTHRESH || res2 > STIFFTHRESH || res3 > STIFFTHRESH ||
        res4 > STIFFTHRESH)
        STIFF = 1;

    return STIFF;
}

void pol_rte_rk4_step(double jI, double jQ, double jU, double jV, double rQ,
                      double rU, double rV, double aI, double aQ, double aU,
                      double aV, double alphaf,double alphai, double dl_current, double C,
                      double complex S_A[]) {
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];
    double delta_a = alphaf-alphai;
    // RK4 with constant coefficients
    // k1
    double complex Ik1 =
        dl_current * C * jI -
        dl_current * C * (aI * I0 + aQ * Q0 + aU * U0 + aV * V0);
    double complex Qk1 =
        dl_current * C * jQ -
        dl_current * C * (aQ * I0 + aI * Q0 + rV * U0 - rU * V0);
    double complex Uk1 =
        dl_current * C * jU -
        dl_current * C * (aU * I0 - rV * Q0 + aI * U0 + rQ * V0);
    double complex Vk1 =
        dl_current * C * jV -
        dl_current * C * (aV * I0 + rU * Q0 - rQ * U0 + aI * V0);

    // k2
    double complex Ik2 = dl_current * C * jI -
                         dl_current * C *
                             (aI * (I0 + 0.5 * Ik1) + aQ * (Q0 + 0.5 * Qk1) +
                              aU * (U0 + 0.5 * Uk1) + aV * (V0 + 0.5 * Vk1));
    double complex Qk2 = dl_current * C * jQ -
                         dl_current * C *
                             (aQ * (I0 + 0.5 * Ik1) + aI * (Q0 + 0.5 * Qk1) +
                              rV * (U0 + 0.5 * Uk1) - rU * (V0 + 0.5 * Vk1));
    double complex Uk2 = dl_current * C * jU -
                         dl_current * C *
                             (aU * (I0 + 0.5 * Ik1) - rV * (Q0 + 0.5 * Qk1) +
                              aI * (U0 + 0.5 * Uk1) + rQ * (V0 + 0.5 * Vk1));
    double complex Vk2 = dl_current * C * jV -
                         dl_current * C *
                             (aV * (I0 + 0.5 * Ik1) + rU * (Q0 + 0.5 * Qk1) -
                              rQ * (U0 + 0.5 * Uk1) + aI * (V0 + 0.5 * Vk1));

    // k3
    double complex Ik3 = dl_current * C * jI -
                         dl_current * C *
                             (aI * (I0 + 0.5 * Ik2) + aQ * (Q0 + 0.5 * Qk2) +
                              aU * (U0 + 0.5 * Uk2) + aV * (V0 + 0.5 * Vk2));
    double complex Qk3 = dl_current * C * jQ -
                         dl_current * C *
                             (aQ * (I0 + 0.5 * Ik2) + aI * (Q0 + 0.5 * Qk2) +
                              rV * (U0 + 0.5 * Uk2) - rU * (V0 + 0.5 * Vk2));
    double complex Uk3 = dl_current * C * jU -
                         dl_current * C *
                             (aU * (I0 + 0.5 * Ik2) - rV * (Q0 + 0.5 * Qk2) +
                              aI * (U0 + 0.5 * Uk2) + rQ * (V0 + 0.5 * Vk2));
    double complex Vk3 = dl_current * C * jV -
                         dl_current * C *
                             (aV * (I0 + 0.5 * Ik2) + rU * (Q0 + 0.5 * Qk2) -
                              rQ * (U0 + 0.5 * Uk2) + aI * (V0 + 0.5 * Vk2));

    // k4
    double complex Ik4 =
        dl_current * C * jI - dl_current * C *
                                  (aI * (I0 + Ik3) + aQ * (Q0 + Qk3) +
                                   aU * (U0 + Uk3) + aV * (V0 + Vk3));
    double complex Qk4 =
        dl_current * C * jQ - dl_current * C *
                                  (aQ * (I0 + Ik3) + aI * (Q0 + Qk3) +
                                   rV * (U0 + Uk3) - rU * (V0 + Vk3));
    double complex Uk4 =
        dl_current * C * jU - dl_current * C *
                                  (aU * (I0 + Ik3) - rV * (Q0 + Qk3) +
                                   aI * (U0 + Uk3) + rQ * (V0 + Vk3));
    double complex Vk4 =
        dl_current * C * jV - dl_current * C *
                                  (aV * (I0 + Ik3) + rU * (Q0 + Qk3) -
                                   rQ * (U0 + Uk3) + aI * (V0 + Vk3));

    double complex x1 = I0 + 1. / 6. * (Ik1 + 2. * Ik2 + 2. * Ik3 + Ik4);
    double complex x2 = Q0 + 1. / 6. * (Qk1 + 2. * Qk2 + 2. * Qk3 + Qk4);
    double complex x3 = U0 + 1. / 6. * (Uk1 + 2. * Uk2 + 2. * Uk3 + Uk4);
    double complex x4= V0 + 1. / 6. * (Vk1 + 2. * Vk2 + 2. * Vk3 + Vk4);

    // if (fabs(alphaf) >0.5E-7){
    S_A[0] = x1;
    S_A[1] = x2*cos(2*delta_a) + x3*sin(2*delta_a);
    S_A[2] = -x2*sin(2*delta_a) + x3*cos(2*delta_a);
    S_A[3] = x4;
    // }
    // else{
    // S_A[0] = x1;
    // S_A[1] = x2;
    // S_A[2] = x3;
    // S_A[3] = x4;
    // }

}

void pol_rte_trapezoid_step(double jI, double jQ, double jU, double jV,
                            double rQ, double rU, double rV, double aI,
                            double aQ, double aU, double aV,double alphaf,double alphai, double dl_current,
                            double C, double complex S_A[]){
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];

    double delta_a = alphaf-alphai;

    double u11 = 1. + 0.5 * dl_current * C * aI;
    double u12 = 0.5 * dl_current * C * aQ;
    double u14 = 0.5 * dl_current * C * aV;
    double l21 = 0.5 * dl_current * C * aQ / u11;
    double u22 = 1. + 0.5 * dl_current * C * aI - l21 * u12;
    double u23 = 0.5 * dl_current * C * rV;
    double u24 = -l21 * u14;
    double l32 = -0.5 * dl_current * C * rV / u22;
    double u33 = 1. + 0.5 * dl_current * C * aI - l32 * u23;
    double u34 = 0.5 * dl_current * C * rQ - l32 * u24;
    double l41 = 0.5 * dl_current * C * aV / u11;
    double l42 = -l41 * u12 / u22;
    double l43 = (-0.5 * dl_current * C * rQ - l42 * u23) / u33;
    double u44 =
        1. + 0.5 * dl_current * C * aI - l41 * u14 - l42 * u24 - l43 * u34;

    // Construct b-vector.
    double b1 =
        I0 + dl_current * C / 2. * (2. * jI - (aI * I0 + aQ * Q0 + aV * V0));
    double b2 =
        Q0 + dl_current * C / 2. * (2. * jQ - (aQ * I0 + aI * Q0 + rV * U0));
    double b3 =
        U0 + dl_current * C / 2. * (2. * jU - (-rV * Q0 + aI * U0 + rQ * V0));
    double b4 =
        V0 + dl_current * C / 2. * (2. * jV - (aV * I0 - rQ * U0 + aI * V0));

    // Construct y.
    double y1 = b1;
    double y2 = b2 - l21 * y1;
    double y3 = b3 - l32 * y2;
    double y4 = b4 - l41 * y1 - l42 * y2 - l43 * y3;

    // Construct x.
    double x4 = y4 / u44;
    double x3 = (y3 - u34 * x4) / u33;
    double x2 = (y2 - u23 * x3 - u24 * x4) / u22;
    double x1 = (y1 - u12 * x2 - u14 * x4) / u11;

    // if (fabs(alphaf) >0.5E-7){
    S_A[0] = x1;
    S_A[1] = x2*cos(2*delta_a) + x3*sin(2*delta_a);
    S_A[2] = -x2*sin(2*delta_a) + x3*cos(2*delta_a);
    S_A[3] = x4;
    // }
    // else{
    // S_A[0] = x1;
    // S_A[1] = x2;
    // S_A[2] = x3;
    // S_A[3] = x4;
    // }

}


int main(void){
 double complex S[4]={2,1,0,0};
 pol_rte_rk4_step(0,0,0,0,0,0,0,0,0,0,0,0.2,0,1,1,S);
 printf("rotation_delta=+.2: I=%.12g Q=%.12g U=%.12g V=%.12g chi=%.12g\n",creal(S[0]),creal(S[1]),creal(S[2]),creal(S[3]),.5*carg(S[1]+I*S[2]));
 printf("stiff_zero=%d stiff_pure_rV=%d stiff_absorb=%d\n",check_stiffness(0,0,0,0,0,0,0,0,0,0,0,1),check_stiffness(0,0,0,0,0,0,1,0,0,0,0,1),check_stiffness(0,0,0,0,0,0,0,1,0,0,0,1));
 for(int n=1;n<=128;n*=2){
  double complex T[4]={0,0,0,0};
  for(int i=0;i<n;i++)pol_rte_rk4_step(1,1,0,0,0,0,0,0,0,0,0,(double)(i+1)/n,(double)i/n,1./n,1,T);
  double complex exact=(1.-cexp(-2.*I))/(2.*I);
  printf("source_split n=%d error=%.12g\n",n,cabs(T[1]+I*T[2]-exact));
 }
 // Synthetic forward traversal: current p=3->1, field A_p=.1*p, camera A0=0.
 double A[5]={0,.1,.2,.3,.4},right=0,actual=0;
 for(int p=3;p>0;p--){right+=A[p-1]-A[p];actual+=A[p+1]-A[p];}
 printf("synthetic telescoping correct_delta=%.12g actual_delta=%.12g, chi_code=%.12g\n",right,actual,-actual);
 return 0;
}
