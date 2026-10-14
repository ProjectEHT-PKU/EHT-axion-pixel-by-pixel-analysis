#include "Axion_Radial.h"
int main(void) {
 double a=.99;
 for(int i=0;i<2;i++) {
  double mu=.3+.01*i,maxall=0,maxoutside=0,at5=0,at1p1=0,rplus=1+sqrt(1-a*a);
  int indexall=-1,indexoutside=-1;
  for(int j=0;j<2400;j++) {
   double r=1+.1*j,X[4]={0,r,0,0},value=Function_R(X,a,mu);
   if(!isfinite(value)){printf("NONFINITE %g %g\n",mu,r);continue;}
   if(value>maxall){maxall=value;indexall=j;}
   if(r>=rplus&&value>maxoutside){maxoutside=value;indexoutside=j;}
   if(j==40)at5=value;
   if(j==1)at1p1=value;
  }
  printf("NORMALIZE a=%g mu=%g rplus=%.17g maxall=%.17g at_r=%.17g maxoutside=%.17g at_r=%.17g finalmax=%.17g R5_old=%.17g R5_exterior=%.17g\n",a,mu,rplus,maxall,1+.1*indexall,maxoutside,1+.1*indexoutside,maxoutside/maxall,at5/maxall,at5/maxoutside);
 }
 return 0;
}
