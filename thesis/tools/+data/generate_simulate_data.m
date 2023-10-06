function [data_raw] = generate_simulate_data(angle,parameters)

% Specify the number of random raw response magnitude curves to generate
% The response magnitude curves are in ideal condition, which means without
% noise disturbance.
%
%   Inputs: ANGLE - Measurement angle of simulated data
%           PARAMETERS - a structure contains information of parameters of tuning curve
%
%   Output: DATA_RAW - store simulate data points
if size(angle,1) < size(angle,2)
    angle = angle'; 
end
data_raw.angle = angle;
data_raw.responses = zeros(length(angle),1);

data_raw.responses = parameters.response + parameters.rpref .*exp(-0.5.*angdiff(angle-parameters.theta_pref).^2/parameters.sigma.^2) + parameters.alpha .* parameters.rpref .* exp(-0.5.*angdiff(angle-parameters.theta_pref+180).^2/parameters.sigma .^2);

if size(data_raw,1) < size(data_raw,2)
    data_raw = data_raw'; 
end