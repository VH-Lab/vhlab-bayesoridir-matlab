#include "filehelp.h"
#include "arrayhelp.h"
#include <stdio.h>
#include <stdlib.h>


int main() {

	double *y1_1 = linspace(1,10,10);
	double y1_2[10];
	FILE *stream;

	printf("Should be 1..10 in steps of 1.\n");
	printdoublearray(y1_1,10);
	printf("\n");
	printf("\n");

	/* write a file */

	stream = fopen("filehelp_test_output.bin","wb");
	if (stream==NULL) {
		printf("File error (opening for writing).\n");
		exit(EXIT_FAILURE);
	}
	lendian_fwrite(y1_1, sizeof(double), 10, stream);
	fclose(stream);

	stream = fopen("filehelp_test_output.bin","rb");
	if (stream==NULL) {
		printf("File error (opening for reading).\n");
		exit(EXIT_FAILURE);
	}

	lendian_fread(y1_2,sizeof(double), 10, stream);
	fclose(stream);

	printf("Should be 1..10 in steps of 1.\n");
	printdoublearray(y1_2,10);
	printf("\n");
	printf("\n");

}

