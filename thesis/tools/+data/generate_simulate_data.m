function [data_raw] = generate_simulate_data(angle,parameters)

% Specify the number of random raw response magnitude curves to generate
% The response magnitude curves are in ideal condition, which means without
% noise disturbance.
%
%   Inputs: ANGLE - Measurement angle of simulated data
%           PARAMETERS - a structure contains information of parameters of tuning curve
%
%   Output: DATA_RAW - store simulate data points

data_raw = zeros(length(angle),1);

data_raw = parameters.response + parameters.rpref .*exp(-0.5.*angdiff(angle-parameters.theta_pref).^2/parameters.sigma.^2) + parameters.alpha .* parameters.rpref .* exp(-0.5.*angdiff(angle-parameters.theta_pref+180).^2/parameters.sigma .^2);