//
//  Axion_Radial.h
//  Radial_function
//
//

#ifndef Axion_Radial_h
#define Axion_Radial_h

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#define CIRCLE_TIMES 1000   // Set the infinite series summation accuracy
#define FLAG_L_M 11         // Set the flag of azimuthal number and angular eigenvalue
#define AZIMUTHAL_M 1.0     // Set the azimuthal number of the perturbation
#define LAMBDA_L_M 2.0      // set the angular eigenvalue Lambda_l_m
#define CIRCULATION_A_0 1.0 // Set the first circulation coefficient
#define NUMBER_OF_E 2.718281828459045

struct Inter_variable
{
    double complex r_plus;
    double complex r_minus;
    double complex omega_c;
    double complex sigma;
    double complex q;
    double complex chi;
    double complex c0;
    double complex c1;
    double complex c2;
    double complex c3;
    double complex c4;
    double complex power1;
    double complex power2;
};

double get_max(double R[200]);
double complex get_omega(double mu, double a);
struct Inter_variable get_Intermediate_variables(double M, double a, double complex omega, double mu);
void change_coefficient(int n, double complex *a, double complex *a_previous, struct Inter_variable A);
double get_radius(double X[4]);
double complex Circulation_R(double r, struct Inter_variable A);
double Function_R(double X[4], double a, double mu);
void print_cable(double a, double Mu_min, double Mu_max);
void print_all_cable();
#endif /* Axion_Radial_h */
