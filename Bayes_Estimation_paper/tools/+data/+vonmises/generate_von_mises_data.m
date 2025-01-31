function [data_raw] = generate_von_mises_data(angle,parameters)

% Specify the number of random raw response curves to generate
% The response curves are in ideal condition, which means without
% noise disturbance.
%
%   Inputs: ANGLE - Measurement angle of simulated data
%           PARAMETERS - a structure contains information of parameters of tuning curve
%
%   Output: DATA_RAW - store simulate data points using von mises model
angle = angle(:);
data_raw.angle = angle;
data_raw.responses = zeros(length(angle),1);

data_raw.responses = parameters.rpref .*exp(parameters.sigma .* ( cosd( 2*angdiff(angle-parameters.theta_pref) ) -1 ) );

if size(data_raw.responses,1) < size(data_raw.responses,2)
    data_raw.responses = data_raw.responses(:); 
end