#include <stdio.h>
int main(int argc,char**argv){FILE*f=fopen(argv[1],"r");char a[100],b[100];double d;for(int i=0;i<17;i++){int n=fscanf(f,"%99s %99s %lf",a,b,&d);printf("%d rc=%d %s %s %.10g\n",i,n,a,b,d);if(n!=3)break;}fclose(f);}
