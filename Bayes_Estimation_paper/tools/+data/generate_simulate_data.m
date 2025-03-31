function [data_raw] = generate_simulate_data(angles,parameters)

% Specify the number of random raw response magnitude curves to generate
% The response magnitude curves are in ideal condition, which means without
% noise disturbance.
%
%   Inputs: ANGLE - Measurement angle of simulated data
%           PARAMETERS - a structure contains information of parameters of tuning curve
%
%   Output: DATA_RAW - store simulate data points
angles = angles(:); 
data_raw.angles = angles;
data_raw.responses = zeros(length(angles),1);

data_raw.responses = parameters.response + parameters.rpref .*exp(-0.5.*angdiff(angles-parameters.theta_pref).^2/parameters.sigma.^2) + parameters.alpha .* parameters.rpref .* exp(-0.5.*angdiff(angles-parameters.theta_pref+180).^2/parameters.sigma .^2);

if size(data_raw.responses,1) < size(data_raw.responses,2)
    data_raw.responses = data_raw.responses(:); 
end