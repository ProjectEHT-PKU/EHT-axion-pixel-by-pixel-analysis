//
//  main.c
//  Radial_function
//
//

#include <stdio.h>
#include "Axion_Radial.h"

int main(int argc, const char *argv[])
{

  // double X[4]={1.0,9.64472,3.0,4.0};             //------------------------
  // Set black hole parameters,and coordinates
  double a = 0.9375; // just for test run
  double Mu_min = 0.1;
  double Mu_max = 0.316;
  // double mu = 0.4;                               //-------------------------
  print_cable(a, Mu_min, Mu_max); // This function print a result.txt in the folder result which contains the Interpolation Funtion R_mu(r).
  //  print_all_cable();              //This function print all the possible Interpolation Function R_mu(r) for the l,m in Axion_Radial.h
  // double R_r_abs = Function_R(X, a, mu);
  return 0;
}
