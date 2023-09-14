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

subplot(3,2,1);
plot(output_struct.marginal_likelihood.Rp.values,output_struct.marginal_likelihood.Rp.likelihoods),
xlabel("Rpref(Hz)"),
ylabel("probability of Rpref");
box off;

subplot(3,2,2);
plot(output_struct.marginal_likelihood.Rn.values,output_struct.marginal_likelihood.Rn.likelihoods),
xlabel("alpha"),
ylabel("probability of alpha");
box off;

subplot(3,2,3);
plot(output_struct.marginal_likelihood.theta_pref.values,output_struct.marginal_likelihood.theta_pref.likelihoods),
xlabel("θpref"),
ylabel("probability of θpref");
box off;

subplot(3,2,4);
plot(output_struct.marginal_likelihood.sigma.values,output_struct.marginal_likelihood.sigma.likelihoods),
xlabel("σ"),
ylabel("probability of σ");
box off;

subplot(3,2,5);
plot(output_struct.marginal_likelihood.Rsp.values,output_struct.marginal_likelihood.Rsp.likelihoods),
xlabel("Rsp"),
ylabel("probability of Rsp");
box off;

subplot(3,2,6);
plot(output_struct.descriptors.oi.values,output_struct.descriptors.oi.likelihoods),
xlabel("Orientation Index"),
ylabel("probability of Orientation Index");

box off;