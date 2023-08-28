#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "arrayhelp.h"

/* linspace - generates a linearly spaced array of numbers between two limits

   Inputs:
	x0 (double)    - the number to start with
	x1 (double)    - the number to end with
	n (int)        - the number of points to make
   Output:
	(double *)     - a pointer to an array with the values (length: n)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Example: 
*/


double * linspace(double x0, double x1, int n) {
	int i;
	double *linspacearray = malloc(sizeof(double)*n);
	double steps = 1;

	if (n>1) steps = n-1;

	if (linspacearray != NULL) {
		for (i=0; i<(steps-1); i++) {
			linspacearray[i] = x0 + i*(x1-x0)/((double)steps);
		}
		/* make the last one exactly right, no round-off */
		if (steps>1) linspacearray[n-1] = x1;
	}
	return(linspacearray);
}

/* logspace - generates a logarithmically spaced array of numbers between limits

   Inputs:
	x0    - the number to start with: base^x0
	x1    - the number to end with: base^x1
	n     - the number of points to make
	base  - the base of the logarithm

   Output:
	(double *)     - a pointer to an array with the values (length: n)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Examples:

   y=logspace(-3,2,100,10) // 100 points between 0.001 and 100 
   y=logspace(log10(0.001),log10(100),100,10) // same 
*/ 

double * logspace(double x0, double x1, int n, double base) {
	int i;
	double *logspacearray = malloc(sizeof(double)*n);
	double stepsize = 0;
	double stepvalue = x0;

	if (n>0) stepsize = (x1-x0)/(n-1);

	if (logspacearray != NULL) {
		for (i=0; i<n-1; i++) {
			logspacearray[i] = pow(base,stepvalue);
			stepvalue += stepsize;
		}
		/* make the last one exactly right, no round-off */
		if (n>1) logspacearray[n-1] = pow(base,x1);
	}
	return(logspacearray);
}

/* printdoublearray - print a double array to standard out

   Inputs:
	double printarray[]  - An array of doubles
	int length           - The length of the array

   Output:
	None (void)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Examples:

   y=logspace(-3,2,10,10); // 10 points between 0.001 and 100 
   printdoublearray(y,10); // print the array
   free(y); // free the memory
*/ 


void printdoublearray(double printarray[], int length) {
	int i;

	for (i=0;i<length;i++) {
		printf("%f ",printarray[i]);
	}
}


