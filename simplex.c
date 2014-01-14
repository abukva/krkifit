#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib.h"

double** create_simplex(int n)
{
	int i,j;
	double MAX=1.0; /* in every dimension the lenght vector will be one unit */

	double **simplex=(double**)malloc((n+1)*sizeof(double*)); /* initialize simplex as a 2D matrix */
	for(i=0;i<(n+1);i++)
	{
		simplex[i]=(double*)malloc(n*sizeof(double));
	}

	/* initialize simplex */
	for(i=0;i<(n+1);i++)
	{
		for(j=0;j<n;j++)
		{
			if(i>0 && (j==(i-1)))
				simplex[i][j]=MAX;
			else
				simplex[i][j]=0.0;
		}
	}
	return simplex;
}

double* minimize(double x[], double y[], int n, int k)
{
	double  **simplex=create_simplex(n); /* create a new simplex */
	
	double *function_values=(double*)malloc((n+1)*sizeof(double)); /* function values for a new simplex */
	function_values=call_functions(x,y,simplex,n,k); /* store function values for every vertex of a simplex in an array */
	
	int worst, second_worst, best; /* variable to store indecies for worst, second_word and best vortex */
	double min, max; /* max and min values for function values */

	int i,j; /* indecies for for loop */

	double* centroid=(double*)malloc(n*sizeof(double)); /* centroid */

	double alpha=1, beta=0.5, gamma=2, delta=0.5; /* scalar parameteras */

	double *vertex_r=(double*)malloc(n*sizeof(double)); /* reflection point */
	double f_r; /* function value at the reflection point */

	double *vertex_e=(double*)malloc(n*sizeof(double)); /* expansion point */
	double f_e; /* function value at the expansion point */

	double *vertex_c=(double*)malloc(n*sizeof(double)); /* contraction point */
	double f_c; /* function value at the contraction point */

	int intterations_max=5000;
	int itteration=0;

	double conv_test; /* check every 42 itterations if function has changde less than 0.01% */
	double diff=42.0;

	do
	{
	/* 1. step Ordering*/

	worst=0, second_worst=0, best=0;
	min=function_values[0], max=function_values[0];
	
	for(i=1;i<(n+1);i++)
	{
		if(min>=function_values[i])
		{
			min=function_values[i];
			best=i;
		}
		if(max<=function_values[i])
		{
			max=function_values[i];
			second_worst=worst;
			worst=i;
		}
	}

	/* 2. step Centroid */

	for(i=0;i<n;i++)
		centroid[i]=0.0;

	for (i=0; i<n; i++)
	{
		for(j=0;j<(n+1);j++)
		{
			if(j!=worst)
			{
				centroid[i]+=simplex[j][i];
			}
		}
		centroid[i]=centroid[i]/n;		
	}

	/* 3. Step Transformation */

	/* Reflect */

	for(i=0;i<n;i++)
	{
		vertex_r[i]=centroid[i]+alpha*(centroid[i]-simplex[worst][i]);
	}
	f_r=objective_function(x,y,vertex_r,n,k);

	if(function_values[best]<=f_r && f_r<function_values[second_worst])
	{
		for(i=0;i<n;i++)
		{
			simplex[worst][i]=vertex_r[i];
		}
		function_values[worst]=f_r;
	}
	else if (f_r<function_values[best]) /* Expand */
	{
		for(i=0;i<n;i++)
		{
			vertex_e[i]=centroid[i]+gamma*(vertex_r[i]-centroid[i]);
		}
		f_e=objective_function(x,y,vertex_e,n,k);
		if(f_e<f_r)
		{
			for(i=0;i<n;i++)
			{
				simplex[worst][i]=vertex_e[i];
			}
			function_values[worst]=f_e;
		}
		else if(f_e>=f_r)
		{
			for(i=0;i<n;i++)
			{
				simplex[worst][i]=vertex_e[i];
			}
			function_values[worst]=f_e;
		}
	}
	else if (f_r>=function_values[second_worst]) /* Contract */
	{
		if(function_values[second_worst]<=f_r && f_r<function_values[worst]) /* Outside */
		{
			for(i=0;i<n;i++)
			{
				vertex_c[i]=centroid[i]+beta*(vertex_r[i]-centroid[i]);
			}
			f_c=objective_function(x,y,vertex_c,n,k);
			if(f_c<=f_r)
			{
				for(i=0;i<n;i++)
				{
					simplex[worst][i]=vertex_c[i];
				}
				function_values[worst]=f_c;
			}
		}
		else if(f_r>=function_values[worst]) /* Inside */
		{
			for(i=0;i<n;i++)
			{
				vertex_c[i]=centroid[i]+beta*(simplex[worst][i]-centroid[i]);
			}
			f_c=objective_function(x,y,vertex_c,n,k);
			if(f_c<function_values[worst])
			{
				for(i=0;i<n;i++)
				{
					simplex[worst][i]=vertex_c[i];
				}
				function_values[worst]=f_c;
			}
		}
		else /* Shrink */
		{
			for(i=0;i<(n+1);i++)
			{
				if(i!=best)
				{
					for(j=0;i<n;j++)
					{
						simplex[i][j]=simplex[best][j]+delta*(simplex[i][j]-simplex[best][j]);
					}
				}
				
			}
			function_values=call_functions(x,y,simplex,n,k);
		}
	}	

	if(itteration%42==0)
	{
		if(i>0)
		{
			diff=abs(conv_test-function_values[best]);
		}
		conv_test=function_values[best];
	}

	itteration++;

	} while (itteration<intterations_max || diff>0.0001);
	printf("%d\n",itteration);
	return simplex[best];

}