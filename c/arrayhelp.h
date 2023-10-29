#ifndef ARRAYHELP_H_
#define ARRAYHELP_H_

#include <math.h>
#include <float.h>
#include <stdio.h>

struct named_array {
        char name[255];
        int array_length;
        double* array_values;
};

typedef struct named_array named_array;

double * linspace(double x0, double x1, int n);

double * logspace(double x0, double x1, int n, double base);

void printdoublearray(double printarray[], int length);

named_array named_array_create(char name[], double array_values[], int array_length);

void named_array_delete(named_array *todelete_array);

int named_array_fwrite(FILE *stream, named_array *to_write);

int named_array_fread(FILE *stream, named_array *to_read);

int named_array_h5read(const char *filename, const char *datapath, named_array *to_read);

#endif // ARRAYHELP_H

