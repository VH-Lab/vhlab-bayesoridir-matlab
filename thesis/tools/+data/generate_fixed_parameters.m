function [parameters_structure] = generate_fixed_parameters(Rpref,Rnull,Theta_pref,Sigma,Response)

% Specify the values ​​of each parameter of the tuning curve in direction
% space
%
%   Inputs: RPREF - above-offset response to the preferred direction
%           RNULL - above-offset response to the null direction
%           THETA_PREF - the preferred direction
%           SIGMA - tuning width
%           RESPONSE - constant offset
%
%   Output: PARAMETERS_STRUCTURE - the structure that store each parameter
%           PARAMETERS_STRUCTURE.ALPHA - calculate the coefficient of Rpref and Rnull (Rnull/Rpref)

parameters_structure.rpref = Rpref;
if Rpref == 0
    parameters_structure.alpha = 0;
else
parameters_structure.alpha = Rnull/Rpref;
end
parameters_structure.theta_pref = Theta_pref;
parameters_structure.sigma = Sigma;
parameters_structure.response = Response;
