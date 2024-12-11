function [output_struct,Lik] = grid_proportional_noise(grid_size,data,noise_mdl, options)
% Perform a Bayesian parameter estimation for double gaussian function
% 
% [OUTPUT_STRUCT, LIK] = GRID_PROPORTIONAL_NOISE(GRID_SIZE, DATA, NOISE_MLD, ...)
%
% Given a set of DATA.angles (such as [0:30:360-30]) and a set of responses to each
% angle, where DATA.responses(i) is the response to a grating presented at
% DATA.angles(i), return a Bayesian estimation of model parameters to the function:
%
% R(THETA) = RSP+RP*exp(-angdiff(OP-ANGLES).^2/(2*SIG^2))+RN*exp(-angdiff(180+OP-ANGLES).^2/(2*SIG^2));
%
% Inputs:
%   GRID_SIZE: A structure with fields:
%      Rsp - the values of Rsp 
%      Rp - the values of Rp to examine
%      Alpha - the values of Rn (equal to Rp * Alpha)
%      Sig - the values of sigma'\
%      
%   DATA : a structure with fields:
%      angles - the angles that were used for stimulation
%      responses - the mean responses to each stimulus
%      num_trials - the number of trials of each stimulus
%   NOISE_MLD  : Model of gaussian standard deviation: sigma = (10^(noise_mld(1)) + response^(noise_mld(2)))/sqrt(num_trials)
%
% OUTPUTS:
%   - OUTPUT_STRUCT: a structure with fields:
%       noise_model - details type of fit, offset, and slope
%       other_parameters - independent variable details
%       marginal_likelihods - marginal likelihoods of parameters
%       maximum_likelihoods_parameters - maximum likelihoods of parameters
%       descriptors - descriptors like direction index, circular variance
%   - LIK: the likelihood matrix
%
% This function looks at the available hardware and chooses to use a PC GPU or Apple Metal if available.
% 
% This function also takes optional parameters:
%    'verbose': 0/1 should we print information about our progress? (default 1)
%    'useCPU': 0/1 should we use the CPU even if a GPU is available? (default 0). Much slower.
%
% Example:
%     % this is a research-grade fit; it will take a long time on a CPU if no GPU is available.
%   rsp_values = sort([-logspace(log10(0.1),log10(40),20) logspace(log10(0.1),log10(40),20)]);
%   rp_values = logspace(log10(0.1),log10(150),100);
%   alpha_values = 0:0.05:1;
%   thetap_values = 0:2:359;
%   sig_values = 10:5:90;
%
%   noise_model = [1.3549 1.9182 0.5461];
%
%   angles = [ 0:30:360-30 ];
%   P = [ 0.5 100 0 45 30 ];
%   resps = vis.oridir.doublegaussianfunc(angles,P);
%
% See also:
%    vis.bayes.noise.proportional()
% 

arguments 
	grid_size (1,1) struct
	data (1,1) struct
	noise_mdl double
	options.verbose (1,1) logical = true
	options.useCPU (1,1) logical = false
end

ops = vlt.data.struct2namevaluepair(rmfield(options,'useCPU'));

if ~options.useCPU & ismac & hasMetal(),
	[output_struct,Lik] = vis.bayes.double_gaussian.grid_proportional_noise_metal(grid_size,data,noise_mdl,ops{:});
elseif ~options.useCPU & canUseGPU(),
	[output_struct,Lik] = vis.bayes.double_gaussian.grid_proportional_noise_gpu(grid_size,data,noise_mdl, ops{:});
else,
	[output_struct,Lik] = vis.bayes.double_gaussian.grid_proportional_noise_cpu(grid_size,data,noise_mdl, ops{:});
end;

