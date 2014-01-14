#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

int main(int argc, char const *argv[])
{
	int n=3; /* number of parameters */
	int k=10; /* number of experimental data */
	int i;
	double* result=(double*)malloc(n*sizeof(double));
	double x[10]={1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	double y[10]={2.98, 6.03, 11.05, 17.99, 27.1, 39.05, 50.89, 67.1, 82.77, 103.2};
	result=minimize(x,y,n,k);
	for(i=0;i<n;i++)
	{
		printf("%.4f ",result[i]);
	}
	printf("\n");
	return 0;
}