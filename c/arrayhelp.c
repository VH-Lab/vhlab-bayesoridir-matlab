#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "arrayhelp.h"
#include "filehelp.h"

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
		for (i=0; i<steps; i++) {
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

named_array named_array_create(char name[], double array_values[], int array_length) {
	named_array new_array;
	int i;

	strcpy(new_array.name,name);
	new_array.array_length = array_length;
	new_array.array_values = (double*)malloc(array_length*sizeof(double));
	if (new_array.array_values != NULL) {
		for (i=0;i<array_length;i++) {
			new_array.array_values[i] = array_values[i];
		}
	} else {
		printf("Error making new array.\n");
		exit(EXIT_FAILURE);
	}
	return(new_array);
}

void named_array_delete(named_array *todelete_array) {
	if ((*todelete_array).array_values != NULL) {
		free((*todelete_array).array_values);
		(*todelete_array).array_values = NULL;
		(*todelete_array).array_length = 0;
		strcpy((*todelete_array).name,"\n");
	}
}

int named_array_fwrite(FILE *stream, named_array *to_write) {
	int output;
	size_t countsize;

	output = fprintf(stream, "%s\n%d\n", (*to_write).name, (*to_write).array_length);

	if (output==0) {
		lendian_fwrite(((*to_write).array_values), sizeof(double), (*to_write).array_length, stream);
	}

	return(output);
}


int named_array_fread(FILE *stream, named_array *to_read) {
	char name[255];
	int array_length;
	double *data;
	int output;
	size_t countsize;

	fscanf(stream, "%s\n%d\n", name, &array_length);

	data = malloc(sizeof(double)*array_length);

	lendian_fread(data, sizeof(double), array_length, stream);

	strcpy((*to_read).name,name);
	(*to_read).array_length = array_length;
	(*to_read).array_values = data;

	return(0);
}


