function [parameters_structure] = generate_von_mises_fixed_parameters(Rpref,Theta_pref,Sigma)

% Specify the values ​​of each parameter of the tuning curve in direction
% space
%
%   Inputs: RPREF - response to the preferred direction
%           THETA_PREF - the preferred direction
%           SIGMA - tuning width
%
%   Output: PARAMETERS_STRUCTURE - the structure that store each parameter

parameters_structure.rpref = Rpref;
parameters_structure.theta_pref = Theta_pref;
parameters_structure.sigma = Sigma;
