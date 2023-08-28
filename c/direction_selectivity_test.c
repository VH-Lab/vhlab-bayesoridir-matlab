#include "direction_selectivity.h"
#include <stdio.h>
#include "arrayhelp.h"



int main() {

	double r, cv, dcv, oi, di;

	double angle_array[8] = {0.0,45.0,90.0,135.0,180.0,225.0,270.0,315.0};
	double r_array[8] = {0,0,0,0,0,0,0,0};

	printf("Angdiff(-90) should be 90: %f\n", angdiff(-90.0));
	printf("Angdiff(90) should be 90: %f\n", angdiff(-90.0));
	printf("Angdiff(10.0-359.0) should be 11: %f\n", angdiff(10.0-359.0));
	printf("Angdiff(10.0-360.0) should be 10: %f\n", angdiff(10.0-360.0));


	r = double_gaussian_angle(0,0,10,0,0,30);

	printf("double_gaussian_angle(0,0,10,0,0,30) should be 10: %f\n", r);

	double_gaussian_curve(r_array,angle_array,8,0,10,0,0,30);

	printdoublearray(r_array, 8);

	printf("\n");

	cv = circular_variance(angle_array, r_array, 8);
	dcv = direction_circular_variance(angle_array, r_array, 8);
	oi = orientation_index_double_gaussian(0,10,0,0,30);
	di = direction_index_double_gaussian(0,10,0,0,30);

	printf("CV is %f and DCV is %f, OI is %f, DI is %f\n", cv, dcv, oi, di);
}
