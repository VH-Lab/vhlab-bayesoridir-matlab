function parameter_struct = parameter_space(description)
% PARAMETER_SPACE - create a parameter space stucture for vis.bayes.double_gaussian.grid_proportional_noise
%
% PARAMETER_STRUCT = PARAMETER_SPACE(DESCRIPTION)
%
% Builds a parameter space structure for vis.bayes.double_gaussian.grid_proportional_noise.
%
% The structure has the following fields:
%
% param_struct = struct('Rsp',rsp_values,...
%    'Rp',rp_values,...
%    'Alpha',alpha_values,...
%    'Op',thetap_values,...
%    'Sig',sig_values);
%
% If DESCRIPTION is 'research-grade', then
%   rsp_values = sort([logspace(log10(0.1),log10(40),40)]);
%   rp_values = logspace(log10(0.1),log10(150),100);
%   alpha_values = 0:0.05:1;
%   thetap_values = 0:2:359;
%   sig_values = 10:5:90;
%
% If DESCRIPTION is 'explore' then
%   rsp_values = sort([logspace(log10(0.1),log10(40),40)]);
%   rp_values = logspace(log10(0.1),log10(150),20);
%   alpha_values = 0:0.05:1;
%   thetap_values = 0:5:359;
%   sig_values = 10:5:90;
%

arguments
    description (1,:) char {mustBeMember(description,{'research-grade','explore'})} = 'research-grade';
end

if strcmp(description,'research-grade')
    rsp_values = sort([-logspace(log10(0.1),log10(40),20) logspace(log10(0.1),log10(40),20)]);
    rp_values = logspace(log10(0.1),log10(150),100);
    alpha_values = 0:0.05:1;
    thetap_values = 0:2:359;
    sig_values = 10:5:90;
end

if strcmp(description,'explore')
    rsp_values = sort([-logspace(log10(0.1),log10(40),20) logspace(log10(0.1),log10(40),20)]);
    rp_values = logspace(log10(0.1),log10(150),20);
    alpha_values = 0:0.05:1;
    thetap_values = 0:5:359;
    sig_values = 10:5:90;
end

parameter_struct = struct('Rsp',rsp_values,...
    'Rp',rp_values,...
    'Alpha',alpha_values,...
    'Op',thetap_values,...
    'Sig',sig_values);
