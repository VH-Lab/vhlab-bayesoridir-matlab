#include "arrayhelp.h"
#include <stdio.h>


int main() {

	double *y1_1 = linspace(1,10,10);
	double *y1_2 = linspace(0.5,0.8,10);

	double *y2_1 = logspace(-2,1,10,10);
	double *y2_2 = logspace(-4,1,10,10);

	named_array na1;

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

	na1 = named_array_create("linspace 1..10 array", y1_1, 10);

	printf("The array name is %s.\n", na1.name);
	printf("Should be 1..10 in steps of of 1 \n");
	printdoublearray(na1.array_values,na1.array_length);
	printf("\n");
	printf("\n");

	named_array_delete(&na1);

}

