#ifndef ARRAYHELP_H_
#define ARRAYHELP_H_

#include <math.h>
#include <float.h>

double * linspace(double x0, double x1, int n);

double * logspace(double x0, double x1, int n, double base);

void printdoublearray(double printarray[], int length);

#endif // ARRAYHELP_H

