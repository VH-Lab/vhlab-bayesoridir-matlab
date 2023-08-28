#include "arrayhelp.h"
#include <stdio.h>


int main() {

	double *y1_1 = linspace(1,10,10);
	double *y1_2 = linspace(0.5,0.8,10);

	double *y2_1 = logspace(-2,1,10,10);
	double *y2_2 = logspace(-4,1,10,10);
	

	printf("Should be 1..10 in steps of 1.\n");
	printdoublearray(y1_1,10);
	printf("\n");
	printf("\n");

	printf("Should be 0.5, 0.08 in 10 steps \n");
	printdoublearray(y1_2,10);
	printf("\n");
	printf("\n");

	printf("Should be 10^-2 to 10 in 10 log steps \n");
	printdoublearray(y2_1,10);
	printf("\n");
	printf("\n");

	printf("Should be 10^-4 to 10 10 log steps \n");
	printdoublearray(y2_2,10);
	printf("\n");
	printf("\n");
}

