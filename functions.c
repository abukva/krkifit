#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

double model_function(double x, double parameters[])
{
	return parameters[0]*x*x+parameters[1]*x+parameters[2];
}

double objective_function(double x[], double y[], double parameters[], int n, int k)
{
	double sum=0;
	int i;
	for(i=0;i<k;i++)
	{
		sum+=(pow(y[i]-model_function(x[i],parameters),2));
	}
	return sum;
}

double *call_functions(double x[], double y[], double **simplex, int n, int k)
{
	int i;
	double* function_values=(double*)malloc((n+1)*sizeof(double));
	for(i=0;i<(n+1);i++)
	{
		function_values[i]=objective_function(x,y,simplex[i],n,k);
	}
	return function_values;
}