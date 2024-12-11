function plot_results(output_struct)
% PLOT_RESULTS - plot the results of double gaussian Bayesian parameter estimation
% 
% PLOT_RESULTS(OUTPUT_STRUCT)
%
% Plots the results of a double Gaussian Bayesian parameter estimation.
% 
% Creates a new figure where all of the marginal probabilities are plotted,
% the most likely parameter values, and the descriptors (orientation index,
% direction index, circular variance, direction circular variance).
%


 % marginals
figure;

marg = sort(fieldnames(output_struct.marginal_likelihoods));

for i=1:5,
    subplot(4,3,i);
    s = getfield(output_struct.marginal_likelihoods,marg{i});
    plot(s.values,s.likelihoods);
    ylabel('Likelihood');
    xlabel(marg{i});
    box off;
end;

subplot(4,3,6);
angles = 0:359;
P = [ output_struct.maximum_likelihood_parameters.parameters.Rsp; ...
      output_struct.maximum_likelihood_parameters.parameters.Rp; ...
      output_struct.maximum_likelihood_parameters.parameters.Rn; ...
      output_struct.maximum_likelihood_parameters.parameters.theta_pref; ...
      output_struct.maximum_likelihood_parameters.parameters.sigma];

r = vis.oridir.doublegaussianfunc(angles,P);
plot(angles,r,'r-');
hold on;
plot(output_struct.other_parameters.independent_variable_value,output_struct.other_parameters.mean_responses,'ko');
box off;

v = {'oi','di','cv','dir_cv'};

for i=1:4,
    subplot(4,3,6+i);
    s = getfield(output_struct.descriptors,v{i});
    plot(s.values,s.likelihoods);
    ylabel('Likelihood');
    xlabel(v{i},'interp','none');
    box off;
end;


