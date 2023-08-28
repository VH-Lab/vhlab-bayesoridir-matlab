#ifndef DIRECTION_SELECTIVITY_H_
#define DIRECTION_SELECTIVITY_H_

#include <math.h>
#include <float.h>

double angdiff(double A);
double double_gaussian_angle(double angle, double offset, double Rp, double Rn, double theta_pref, double sigma);
void double_gaussian_curve(double *r, double angle[], int angle_n, double offset, double Rp, double Rn, double theta_pref, double sigma);

double orientation_index_double_gaussian(double offset, double Rp, double Rn, double theta_pref, double sigma);
double direction_index_double_gaussian(double offset, double Rp, double Rn, double theta_pref, double sigma);

double circular_variance(double angle[], double r[], int angle_n);
double direction_circular_variance(double angle[], double r[], int angle_n);

double direction_ratio_help(double numerator, double denominator);

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
        double noise_model_offset, double noise_model_sigma);


#endif // DIRECTION_SELECTIVITY_H_

