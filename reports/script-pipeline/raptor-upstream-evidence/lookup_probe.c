#include <stdio.h>
#include <stdlib.h>
int main(int argc,char**argv){double table[23][2400][3]={0};for(int i=0;i<23;i++)table[i][2][0]=.09+.01*i;double a=atof(argv[1]);int i=1;while(table[i][2][0]<a)++i;--i;printf("alpha=%g bracket=%d,%d\n",a,i,i+1);return 0;}
