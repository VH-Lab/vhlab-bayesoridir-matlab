function [output_struct,Lik,debug] = grid_proportional_noise_metal(param_grid, resp_struct, noise_mdl, options)
% GRID_PROPRTIONAL_NOISE_METAL - calculates Bayesian likelihood of double gaussian model using Apple Metal
%
% [OUTPUT_STRUCT, LIK]  = GRID_PROPORTIONAL_NOISE_METAL(PARAM_GRID, RESP_STRUCT, NOISE_MDL)
%
% Given a set of RESP_STRUCT.angles (such as [0:30:360-30]) and a set of responses to each
% angle, where RESP_STRUCT.mean_responses(i) is the mean response to a grating presented at
% RESP_STRUCT.angles(i), return a Bayesian estimation of model parameters to the function:
%
% R(THETA) = RSP+RP*exp(-vlt.math.angdiff(OP-ANGLES).^2/(2*SIG^2))+RN*exp(-vlt.math.angdiff(180+OP-ANGLES).^2/(2*SIG^2));
%
% - Response variance is proportional to mean response for pooled data in log space
%
% INPUTS:
%   - PARAM_GRID: A structure with fields:
%       Rsp - the values of Rsp to examine
%       Rp - the values of Rp to examine
%       Alpha - the values of Rn to examine (Rn = Rp*Alpha)
%       Op - the values of OP to examine 
%       Sig - the values of SIG to examine
%   - RESP_STRUCT: a structure with fields:
%       angle - the contrast that were used for stimulation
%       mean_responses - the mean responses to each stimulus
%       num_trials - the number of trials of each stimulus
%   - NOISE_MDL: a log linear regression model of noise for the given data [offset slope] or [offset k slope]
%
% OUTPUTS:
%   - OUTPUT_STRUCT: a structure with fields:
%       noise_model - details type of fit, offset, and slope
%       other_parameters - independent variable details
%       marginal_likelihods - marginal likelihoods of parameters
%       maximum_likelihoods_parameters - maximum likelihoods of parameters
%       descriptors - descriptors like direction index, circular variance
%   - LIK
%
% The prior probabilities of the values of the parameters are assumed
% to be uniform over the values called 
%
% See also:
%   vis.bayes.noise.proportional()

arguments
        param_grid (1,1) struct
        resp_struct (1,1) struct
        noise_mdl double
        options.verbose (1,1) logical = true
end


%% DEFINE VARIABLES
% from data_struct
resps = resp_struct.mean_responses(:);
angle_value = resp_struct.angles(:);
num_trials = resp_struct.num_trials(:);

% from grid_size
rsp_values = param_grid.Rsp;
rp_values = param_grid.Rp;
alpha_values = param_grid.Alpha;
op_values = param_grid.Op;
sig_values = param_grid.Sig;

% from noise_mdl
if numel(noise_mdl)>2,
	noise_mdl = [2;noise_mdl(:)];
else,
	noise_mdl = [1;noise_mdl(:)];
end;

funcName = 'doubleGaussianProportionalNoiseLog';
here = which('vis.bayes.double_gaussian.gpudemo');
[parent,me] = fileparts(here);
kernel = fileread(fullfile(parent,'doublegaussianproportionalnoise.mtl'));

lik = zeros(numel(rsp_values)*numel(rp_values)*numel(alpha_values)*numel(op_values)*numel(sig_values),1,'single');
histogram_index = zeros(numel(rsp_values)*numel(rp_values)*numel(alpha_values)*numel(op_values)*numel(sig_values),1,'single');

oi_bins  = 0:0.05:1.00;
di_bins  = 0:0.05:1.00;
cv_bins  = 0:0.05:1.00;
dcv_bins = 0:0.05:1.00;

if numel(num_trials)~=numel(angle_value),
	error(['There must be a numtrial entry for each response. The size of RESP_STRUCT.num_trials must match RESP_STRUCT.mean_responses.']);
end;

angles_and_response = single([numel(angle_value) ; angle_value(:); resps(:); num_trials(:)]);
grid = single([numel(rsp_values);
        numel(rp_values);
        numel(alpha_values);
        numel(op_values);
        numel(sig_values);
        rsp_values(:);
        rp_values(:);
        alpha_values(:);
        op_values(:);
        sig_values(:)]);

bin_values = single([numel(oi_bins); numel(di_bins); numel(cv_bins); numel(dcv_bins); ...
  oi_bins(:); di_bins(:); cv_bins(:); dcv_bins(:)]);

answers = single(zeros(100,1));

metalConfig = MetalConfig;
device = MetalDevice(metalConfig.gpudevice);

bufferLik = MetalBuffer(device,lik);
bufferAnglesAndResponse = MetalBuffer(device,angles_and_response);
bufferGrid = MetalBuffer(device,grid);
bufferNoiseModel = MetalBuffer(device,single(noise_mdl));
bufferBinValues = MetalBuffer(device,bin_values);
bufferHistogramIndex = MetalBuffer(device,histogram_index);

bufferAnswers = MetalBuffer(device,answers);

%disp(['About to call Metal'])

MetalCallKernel(funcName,{bufferLik,bufferAnglesAndResponse,bufferGrid,bufferNoiseModel,bufferBinValues,bufferHistogramIndex,bufferAnswers},kernel);

%disp(['Finished with Metal Call'])

debug.answers = single(bufferAnswers);

lik = single(bufferLik);
lik = reshape(lik,numel(rsp_values),numel(rp_values),numel(alpha_values),numel(op_values),numel(sig_values));

histogram_index = uint64(single(bufferHistogramIndex));

debug.histogram_index = histogram_index;

%disp(['Done pulling variable info'])

Lik = 10.^double(lik);
Lik_prenorm = Lik;
Lik = Lik./sum(Lik(:));

%% EXTRACT MARGINAL LIKELIHOOD OF EACH PARAMETER 
lik_rsp = squeeze(sum(sum(sum(sum(Lik,5),4),3),2));
lik_rsp = lik_rsp./sum(lik_rsp,"all");

lik_rp = squeeze(sum(sum(sum(sum(Lik,5),4),3),1))';
lik_rp = lik_rp./sum(lik_rp,"all");

lik_alpha = squeeze(sum(sum(sum(sum(Lik,5),4),2),1));
lik_alpha = lik_alpha./sum(lik_alpha,"all");

lik_op = squeeze(sum(sum(sum(sum(Lik,5),3),2),1));
lik_op = lik_op./sum(lik_op,"all");

lik_sig = squeeze(sum(sum(sum(sum(Lik,4),3),2),1));
lik_sig = lik_sig./sum(lik_sig,"all");

%% MAXIMUM LIKELIHOOD CONDITIONS
[~,ind] = max(Lik,[],'all');
[rsp_ml,rp_ml,alpha_ml,op_ml,sig_ml] = ind2sub(size(Lik),ind);
dgresp = vis.oridir.doublegaussianfunc(angle_value,...
	[rsp_values(rsp_ml) rp_values(rp_ml) alpha_values(alpha_ml) op_values(op_ml) sig_values(sig_ml)]);

%% DESCRIPTOR estimation
%disp('Descriptor calculations')

oi_lik = 0 * oi_bins;
di_lik = 0 * di_bins;
cv_lik = 0 * cv_bins;
dcv_lik = 0 * dcv_bins;

maxIndexes = 1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)];

for i=1:4,
	indexes = vlt.data.ind2subD(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],histogram_index,i);
	for bin = 1:maxIndexes(i)
		if i==1,
			oi_lik(bin) = sum(Lik(indexes==bin));
		elseif i==2,
			di_lik(bin) = sum(Lik(indexes==bin));
		elseif i==3,
			cv_lik(bin) = sum(Lik(indexes==bin));
		elseif i==4,
			dcv_lik(bin) = sum(Lik(indexes==bin));
		end;
	end;
end;

oi_lik = oi_lik./sum(oi_lik);
di_lik = di_lik./sum(di_lik);
cv_lik = cv_lik./sum(cv_lik);
dcv_lik = dcv_lik./sum(dcv_lik);

%% CREATE OUTPUT_STRUCT
output_struct = struct( ...
	'noise_model',struct('type',{'proportional'},'noise_model',noise_mdl), ...
	'other_parameters',struct('independent_variable',{'angle'},'independent_variable_value',resp_struct.angles,...
		'mean_responses',resp_struct.mean_responses,'num_trials',resp_struct.num_trials), ...
	'marginal_likelihoods',struct( ...
		'Rsp',struct('values',param_grid.Rsp,'likelihoods',lik_rsp), ...
		'Rp',struct('values',param_grid.Rp,'likelihoods',lik_rp), ...
		'Alpha',struct('values',param_grid.Alpha,'likelihoods',lik_alpha), ...
		'theta_pref',struct('values',param_grid.Op,'likelihoods',lik_op), ...
		'sigma',struct('values',param_grid.Sig,'likelihoods',lik_sig)), ...
	'maximum_likelihood_parameters',struct(...
			'parameters',...
				struct('Rsp',param_grid.Rsp(rsp_ml),...
				'Rp',param_grid.Rp(rp_ml),...
				'Rn',param_grid.Alpha(alpha_ml)*param_grid.Rp(rp_ml),...
				'theta_pref',param_grid.Op(op_ml),...
				'sigma',param_grid.Sig(sig_ml)),...
			'tuning_curve',dgresp),...
	'descriptors', struct(...
		'oi', struct('values',[oi_bins NaN],'likelihoods',oi_lik), ...
		'di', struct('values',[di_bins NaN],'likelihoods',di_lik), ...
		'cv', struct('values',[cv_bins NaN],'likelihoods',cv_lik), ...
		'dir_cv', struct('values',[dcv_bins NaN],'likelihoods',dcv_lik) ) );

