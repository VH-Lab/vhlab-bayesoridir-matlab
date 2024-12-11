function plot(output_struct)
% PLOT - plot the output of the Bayesian Double Gaussian grid model
%
% vis.bayes.double_gaussian.plot(OUTPUT_STRUCT)
%
% Plot the output of the Bayesian double Gaussian fit in the current
% figure.
%
% See also: vis.bayes.double_gaussian.grid_proportional_noise()
%

subplot(3,3,1);
plot(output_struct.marginal_likelihood.Rp.values,output_struct.marginal_likelihood.Rp.likelihoods,'b-o'),
set(gca,'xscale','log')
xlabel("Rpref(Hz)"),
ylabel("probability of Rpref");
box off;

subplot(3,3,2);
plot(output_struct.marginal_likelihood.Alpha.values,output_struct.marginal_likelihood.Alpha.likelihoods,'b-o'),
xlabel("alpha"),
ylabel("probability of alpha");
box off;

subplot(3,3,3);
plot(output_struct.marginal_likelihood.theta_pref.values,output_struct.marginal_likelihood.theta_pref.likelihoods,'b-o'),
xlabel("θpref"),
ylabel("probability of θpref");
box off;

subplot(3,3,4);
plot(output_struct.marginal_likelihood.sigma.values,output_struct.marginal_likelihood.sigma.likelihoods,'b-o'),
xlabel("σ"),
ylabel("probability of σ");
box off;

subplot(3,3,7);
plot(output_struct.marginal_likelihood.Rsp.values,output_struct.marginal_likelihood.Rsp.likelihoods,'b-o'),
set(gca,'xscale','log')
xlabel("Rsp"),
ylabel("probability of Rsp");
box off;

subplot(3,3,5);
plot(output_struct.descriptors.oi.values,output_struct.descriptors.oi.likelihoods,'b-o'),
xlabel("Orientation Index"),
ylabel("probability of Orientation Index");

subplot(3,3,8);
plot(output_struct.descriptors.di.values,output_struct.descriptors.di.likelihoods,'b-o'),
xlabel("Direction Index"),
ylabel("probability of Direction Index");


subplot(3,3,6);
plot(output_struct.descriptors.cv.values,output_struct.descriptors.cv.likelihoods,'b-o'),
xlabel("CV"),
ylabel("Circular variance (orientation space)");

subplot(3,3,9);
plot(output_struct.descriptors.dir_cv.values,output_struct.descriptors.dir_cv.likelihoods,'b-o'),
xlabel("DCV"),
ylabel("probability of DCV");

box off;