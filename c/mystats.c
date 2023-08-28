#include "mystats.h"
#include <math.h>
#include <float.h>

const double MYSTATS_1_OVER_SQRT_2PI = 0.398942280401433;

/* normpdf - The probability density function of the normal distribution at a
             given value of x.

   Inputs:
	x (double)                  - The value to evaluate
        mean (double)               - the mean of the distribution
	standard_deviation(double)  - the standard deviation 
   Output:
        double                      - the value of the PDF

   Example:
   double p = normpdf(0,0,1); // prob that x==0 for a dist w/ mean 0, std 1
*/

double normpdf(double x, double mean, double standard_deviation) {
	double myvalue;
	return( (MYSTATS_1_OVER_SQRT_2PI/standard_deviation) * (exp(-(x-mean)*(x-mean)/(2*standard_deviation*standard_deviation))) );
}
