#include "direction_selectivity.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "mystats.h"

/* angdiff - Angular difference in 0..360

double D = angdiff(double A)

Computes the angular difference in 0..360. Answer is returned in degrees. 

Returns min(fabs([A;A+360;A-360]));

  Inputs:
	A (double) - the angular difference
  Outputs:
	(double)   - the angular difference wrapped to 0..360

Examples:
  angdiff(0-90) = 90 // there is a 90 degree angular difference between 0, 90
  angdiff(90-0) = 90 // there is a 90 degree angular difference between 90, 0
  angdiff(10-359) = 11  // 11 degree angular difference between 10, 359
  angdiff(10-360) = 10 // 10 degree angular difference between 10 , 360

*/

double angdiff(double A) {
	return(fmin(fmin( fabs(A), fabs(A+360)), fabs(A-360)));
}

/* double_gaussian_angle - double gaussian calculation in angular space

Computes r = offset + Rp*exp(- (angdiff(angle-theta_pref)^2)/(2*sigma^2)) + Rn*exp(- (angdiff(angle-theta_pref-180)^2)/(2*sigma^2)

Inputs: 
	angle - the angle, in degrees
	offset - the offset of the response from 0 that does depend on angle
	Rp - the response above offset at the preferred angle theta_pref
	Rn - the response above offset at the opposite angle theta_pref+180
	theta_pref - the preferred angle, in degrees
	sigma - the tuning width of the double gaussiana

Outputs:
	(double) - r, the response in the equation above

Example:
  r = double_gaussian_angle(0, 0, 10, 0, 0, 20); // should be 10

*/

double double_gaussian_angle(double angle, double offset, double Rp, double Rn, double theta_pref, double sigma) {
	double gdenom = 1.0/(2*sigma*sigma);
	double gpref = angdiff(angle-theta_pref);
	double gnull = angdiff(angle-(theta_pref+180));
	
	return( offset+Rp*exp(-gpref*gpref*gdenom)+Rn*exp(-gnull*gnull*gdenom) );
}


/* double_gaussian_curve - compute the responses to a full double gaussian curve in angular space

Computes the response to a full double gaussian curve with angle values in an array

Inputs:
	r - a pointer to a double array, which must already be allocated and be the same size as the angle array
	angle - an array of angle values
	angle_n - the number of entries in the array of angle values
	offset - the offset of the response from 0 that does depend on angle
	Rp - the response above offset at the preferred angle theta_pref
	Rn - the response above offset at the opposite angle theta_pref+180
	theta_pref - the preferred angle, in degrees
	sigma - the tuning width of the double gaussiana

Outputs:
	none

*/

void double_gaussian_curve(double *r, double angle[], int angle_n, double offset, double Rp, double Rn, double theta_pref, double sigma) {
	int a;

	for (a=0;a<angle_n;a++) {
		r[a] = double_gaussian_angle(angle[a],offset,Rp,Rn,theta_pref,sigma);
	}
}

/* circular_variance - circular variance in orientation space

Computes the circular variance in orienatation space (Mazurek et al. 2014)

If the responses are all 0, then the circular variance is defined to be 1.

Inputs:
	angle - an array of angle values
	r - an array of response values
	angle_n - the number of entries in the array of angle values

Outputs:
	the circular variance

*/

double circular_variance(double angle[], double r[], int angle_n) {
	int a;
	double vec_real = 0.0, vec_complex = 0.0, total = 0.0, angle_rad, numerator, L;

	for (a=0;a<angle_n;a++) {
		angle_rad = 2*remainder(angle[a],180.0)*M_PI/180;
		vec_real += r[a]*cos(angle_rad);
		vec_complex += r[a]*sin(angle_rad);
		total += r[a];
	}
	if (total == 0) {
		return(1.0);
	}
	numerator = sqrt(vec_real*vec_real+vec_complex*vec_complex);
	L = numerator/total;
	return(1.0 - L);
}

/* direction_circular_variance - circular variance in direction space

Computes the circular variance in direction space (Mazurek et al. 2014)

If the responses are all 0, the direction circular variance is 1.0.

Inputs:
	angle - an array of angle values
	r - an array of response values
	angle_n - the number of entries in the array of angle values

Outputs:
	the circular variance in direction space

*/

double direction_circular_variance(double angle[], double r[], int angle_n) {
	int a;
	double vec_real = 0, vec_complex = 0, total = 0, angle_rad, numerator, L;

	for (a=0;a<angle_n;a++) {
		angle_rad = angle[a]*M_PI/180;
		vec_real += r[a]*cos(angle_rad);
		vec_complex += r[a]*sin(angle_rad);
		total += r[a];
	}
	if (total==0) {
		return(1.0);
	}
	numerator = sqrt(vec_real*vec_real+vec_complex*vec_complex);
	L = numerator/total;
	return(1.0 - L);
}

/* orientation index - orientation index

Computes the orientation index for a double gaussian fit of orientation/direction responses (Mazurek et al. 2014)

OI =  (R(theta_pref)+R(theta_pref+180)-R(theta_pref+90)-R(theta_pref-90))/(R(theta_pref)+R(theta_pref+180)

If Rpref+Rnull is 0, the OI index 0.

Inputs:
	offset - the offset of the response from 0 that does depend on angle
	Rp - the response above offset at the preferred angle theta_pref
	Rn - the response above offset at the opposite angle theta_pref+180
	theta_pref - the preferred angle, in degrees
	sigma - the tuning width of the double gaussiana

Outputs:
	the circular variance in direction space

*/

double orientation_index_double_gaussian(double offset, double Rp, double Rn, double theta_pref, double sigma) {
	double rpref = double_gaussian_angle(theta_pref,offset,Rp,Rn,theta_pref,sigma);
	double rnull = double_gaussian_angle(theta_pref+180,offset,Rp,Rn,theta_pref,sigma);
	double rorth_plus = double_gaussian_angle(theta_pref+90,offset,Rp,Rn,theta_pref,sigma);
	double rorth_minus = double_gaussian_angle(theta_pref-90,offset,Rp,Rn,theta_pref,sigma);

	if ( (rpref+rnull) == 0) return (0);
	return( (rpref+rnull-rorth_plus-rorth_minus)/(rpref+rnull));
}

/* direction index - direction index

Computes the direction index for a double gaussian fit of orientation/direction responses (Mazurek et al. 2014)

DI =  (R(theta_pref)-R(theta_pref+180))/R(theta_pref)

If Rpref is 0, the DI index 0.

Inputs:
	offset - the offset of the response from 0 that does depend on angle
	Rp - the response above offset at the preferred angle theta_pref
	Rn - the response above offset at the opposite angle theta_pref+180
	theta_pref - the preferred angle, in degrees
	sigma - the tuning width of the double gaussiana

Outputs:
	the direction index DI

*/

double direction_index_double_gaussian(double offset, double Rp, double Rn, double theta_pref, double sigma) {
	double rpref = double_gaussian_angle(theta_pref,offset,Rp,Rn,theta_pref,sigma);
	double rnull = double_gaussian_angle(theta_pref+180,offset,Rp,Rn,theta_pref,sigma);
	
	if (rpref == 0) return (0);

	return((rpref-rnull)/rpref);
}


/* direction_ratio_help - return the ratio of two numbers, or NAN, INFINITY, or -INFINITY as appropriate

If the denominator is 0, then NAN is returned if the numerator is also 0, or INFINITY if the numerator is positive
and -INFINITY is the numerator is negative.

Inputs:
	numerator - the numerator
	denominator - the denominator
Outputs:
	the ratio, adjusted for not-a-number or infinity as needed

*/

double direction_ratio_help(double numerator, double denominator) {
	if (denominator==0) {
		if (numerator==0) {
			return(NAN);
		} else {
			if (numerator>0) {
				return(INFINITY);
			} else {
				return(-INFINITY);
			}
		}
	} else {
		return(numerator/denominator);
	}
}


/* double_gaussian_bayes - compute marginal likelihoods of double gaussian parameters


 // note: right now does not compute oi, di, cv, dcv marginals


Inputs:
	angle_values_num - the number of angles tested
	angle_values - the angle values tested
	response_values - the mean response at each angle tested,
	number_of_measurements - the number of measurements at each angle tested
	offset_values_num - the number of elements in the double gaussian offset array
	offset_values - an array of double gaussian fit offset values to examine
	offset_values_prior - the prior probability of each offset value
	offset_marginal - the marginal probability of each offset_value (returned back), should already be created
	Rp_values_num - the number of elements in the double gaussian Rp array
	Rp_values - an array of double gaussian fit Rp values to examine
	Rp_values_prior - the prior probability of each Rp value
	Rp_marginal - the marginal probability of each Rp_value (returned back), should already be created
	Rn_values_num - the number of elements in the double gaussian Rn array
	Rn_values - an array of double gaussian fit Rn values to examine
	Rn_values_prior - the prior probability of each Rn value
	Rn_marginal - the marginal probability of each Rn_value (returned back), should already be created
	theta_pref_values_num - the number of elements in the double gaussian theta_pref array
	theta_pref_values - an array of double gaussian fit theta_pref values to examine
	theta_pref_values_prior - the prior probability of each theta_pref value
	theta_pref_marginal - the marginal probability of each theta_pref_value (returned back), should already be created
	sigma_values_num - the number of elements in the double gaussian sigma array
	sigma_values - an array of double gaussian fit sigma values to examine
	sigma_values_prior - the prior probability of each sigma value
	sigma_marginal - the marginal probability of each sigma_value (returned back), should already be created
	oi_bin_number - the number of bin edges for oi
	oi_bin_edges - binned values over which to compute oi probability
	oi_marginal - the marginal probabily of seeing oi for each value, should already be created
	di_bin_number - the number of bin edges for di
	di_bin_edges - binned values over which to compute di probability
	di_marginal - the marginal probabily of seeing di for each value, should already be created
	cv_bin_number - the number of bin edges for cv
	cv_bin_edges - binned values over which to compute cv probability
	cv_marginal - the marginal probabily of seeing dcv for each value, should already be created
	dcv_bin_number - the number of bin edges for dcv
	dcv_bin_edges - binned values over which to compute dcv probability
	dcv_marginal - the marginal probabily of seeing dcv for each value, should already be created
	noise_model_offset - the noise model offset
	noise_model_sigma - the noise model sigma
	
Output:
	returns 0 if successful, 1 if memory could not be created for likelihood functions
*/ 

int double_gaussian_bayes(
	int angle_values_num, double angle_values[], double response_values[], int number_of_measurements[],
	int offset_values_num, double offset_values[], double offset_values_prior[], double *offset_marginal,
	int Rp_values_num, double Rp_values[], double Rp_values_prior[], double *Rp_marginal,
	int Rn_values_num, double Rn_values[], double Rn_values_prior[], double *Rn_marginal,
	int theta_pref_values_num, double theta_pref_values[], double theta_pref_values_prior[], double *theta_pref_marginal,
	int sigma_values_num, double sigma_values[], double sigma_values_prior[], double *sigma_marginal,
	int oi_bin_number, double oi_bin_edges[], double *oi_marginal,
	int di_bin_number, double di_bin_edges[], double *di_marginal,
	int cv_bin_number, double cv_bin_edges[], double *cv_marginal,
	int dcv_bin_number, double dcv_bin_edges[], double *dcv_marginal,
	double noise_model_offset, double noise_model_sigma) { 

	long int offset_, Rp_, Rn_, theta_pref_, sigma_, a_, i, lik_N, counter = 0;
	long int offset_step, Rp_step, Rn_step, theta_pref_step, sigma_step;

	double *resp_here, *lik, prob_here, lik_sum;

	resp_here = (double *) malloc(sizeof(double)*angle_values_num);
	if (resp_here == NULL) return(1); 

	sigma_step = 1;
	theta_pref_step = sigma_values_num;
	Rn_step = theta_pref_values_num*theta_pref_step;
	Rp_step = Rn_values_num*Rn_step;
	offset_step = Rp_values_num*Rp_step;

	lik_N = offset_values_num*Rp_values_num*Rn_values_num*theta_pref_values_num*sigma_values_num;

	lik = (double *)malloc(sizeof(double)*lik_N);
	if (lik == NULL) return(1); // failed to allocate memory

	for (offset_=0; offset_<offset_values_num; offset_++) {
	  for (Rp_=0;Rp_<Rp_values_num;Rp_++) {
	    for (Rn_=0;Rn_<Rn_values_num;Rn_++) {
	      for (theta_pref_=0;theta_pref_<theta_pref_values_num;theta_pref_++) {
	        for (sigma_=0;sigma_<sigma_values_num;sigma_++) {
			double_gaussian_curve(resp_here, angle_values, angle_values_num, 
				offset_values[offset_], Rp_values[Rp_], Rn_values[Rn_], theta_pref_values[theta_pref_], sigma_values[sigma_]);
			prob_here = 1.0;
			for (a_=0;a_<angle_values_num;a_++) {
				prob_here *= normpdf(resp_here[a_]-response_values[a_], 0.0, 
					noise_model_offset+noise_model_sigma*resp_here[a_]/sqrt(number_of_measurements[a_]));
			}
			lik[counter] = prob_here;
			counter+=1;
			lik_sum += prob_here;
	        }
	      }
	    }
	  }
	}

	// normalize
	for (i=0;i<lik_N;i++) lik[i] = lik[i]/lik_sum;

	// compute marginals
	for (offset_=0; offset_<offset_values_num; offset_++) {
	  for (Rp_=0;Rp_<Rp_values_num;Rp_++) {
	    for (Rn_=0;Rn_<Rn_values_num;Rn_++) {
	      for (theta_pref_=0;theta_pref_<theta_pref_values_num;theta_pref_++) {
	        for (sigma_=0;sigma_<sigma_values_num;sigma_++) {
			offset_marginal[offset_] += lik[sigma_+theta_pref_*theta_pref_step+Rn_*Rn_step+Rp_*Rp_step+offset_*offset_step];
			Rp_marginal[Rp_] += lik[sigma_+theta_pref_*theta_pref_step+Rn_*Rn_step+Rp_*Rp_step+offset_*offset_step];
			Rn_marginal[Rn_] += lik[sigma_+theta_pref_*theta_pref_step+Rn_*Rn_step+Rp_*Rp_step+offset_*offset_step];
			theta_pref_marginal[theta_pref_] += lik[sigma_+theta_pref_*theta_pref_step+Rn_*Rn_step+Rp_*Rp_step+offset_*offset_step];
			sigma_marginal[sigma_] += lik[sigma_+theta_pref_*theta_pref_step+Rn_*Rn_step+Rp_*Rp_step+offset_*offset_step];
	        }
	      }
	    }
	  }
	}

	free(resp_here);
	free(lik);

	return(0);
}
