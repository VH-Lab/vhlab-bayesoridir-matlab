#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>



int main() {
	double i,j,k;
	double a=1.05;
	double max_value = 1000;

	


	clock_t start, end;
	double cpu_time_used;
     
	printf("This is a test.\nI haven't written one of these in awhile.\n");

	start = clock();
	for (i=0; i<max_value; i++) {
		for (j=0; j<max_value; j++) {
			for (k=0; k<max_value; k++) {
				a = 5*a + 12 - 32 + 55.0340 - 5*a + exp(-a) - exp(-a) + 2*i+3*j+4*k;
			}
		}
	}

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	printf("Operations: %e\n", (double)(max_value*max_value*max_value));
	printf("CPU time used: %f\n", cpu_time_used);

	return(0);
}

/* linspace - generates a linearly spaced set of numbers between two limits

   Inputs:
	n0 (double)    - the number to start with
	n1 (double)    - the number to end with
	steps (int)    - the number of steps to use
   Output:
	(double *)     - a pointer to an array with the values (length: steps)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Example: 
*/


(double*) linspace(double n0, double n1, int steps) {

	int i;
	double *linspacearray = malloc(sizeof(double)*steps);

	if (linspacearray != NULL) {
		for (i=0; i<steps; i++) {
			linspacearray[i] = (n1-n0)/((double)steps);
		}
	}
}

/* logspace - generates a logarithmically spaced set of double numbers between two limits

   Inputs:
	n0    - the number to start with
	n1    - the number to end with
	steps - the number of steps to use
	base  - the base of the logarithm




*/ 

(double *) logspace(double n0, double n1, int steps, double base)
