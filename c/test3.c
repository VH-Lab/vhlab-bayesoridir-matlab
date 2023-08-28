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


