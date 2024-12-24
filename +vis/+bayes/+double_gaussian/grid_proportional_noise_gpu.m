function [output_struct,Lik] = grid_proportional_noise_gpu(grid_size,data,noise_mdl,varargin)
%   DESCRIPTION:
%
%   Inputs: GRID_SIZE - 
%           DATA - 
%           NOISE_MDL - 
%
%   Outputs:OUTPUT_STRUCT - 
%           LIK - 
%
% extract coefficients from noise linear regression model


if numel(noise_mdl) == 2
    offset = noise_mdl(1);
    slope = noise_mdl(2);
    k = [];
elseif numel(noise_mdl) == 3
    offset = noise_mdl(1);
    k = noise_mdl(2);
    slope = noise_mdl(3);
end
% build bayes grid matrix
Lik = gpuArray(zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp)));
ori_cv_value = Lik;
dir_cv_value = Lik;
[rp_expanded,op_expanded,alpha_expanded,sig_expanded] = ndgrid(gpuArray(grid_size.Rp),gpuArray(grid_size.Op),gpuArray(grid_size.Alpha),gpuArray(grid_size.Sig));
data_angle_expanded = gpuArray(reshape(data.angles, 1, 1, 1, 1, []));
data_mean_expanded = gpuArray(reshape(data.mean_responses', 1, 1, 1, 1, []));
num_trials = gpuArray(reshape(data.num_trials', 1, 1, 1, 1, []));
num_trials_expanded = repmat([num_trials],[size(rp_expanded) 1]);

for rsp = 1:length(grid_size.Rsp)
    fprintf('This loop round left: %d \n',length(grid_size.Rsp)+1-rsp);
    vrsp = grid_size.Rsp(rsp) + rp_expanded .* exp(-0.5 * angdiff(data_angle_expanded - op_expanded) .^ 2 ./ sig_expanded .^ 2) + alpha_expanded .* rp_expanded .* exp(-0.5 * angdiff(data_angle_expanded - (op_expanded + 180)) .^ 2 ./ sig_expanded .^ 2);
    noise_sigma = vis.bayes.noise.proportional(noise_mdl,abs(vrsp),num_trials_expanded);
    prsp = normpdf(data_mean_expanded,vrsp,noise_sigma);
    multiprsp = squeeze(prod(prsp,5));
    Lik(:,:,:,:,rsp) = multiprsp;
    
    %circular variance and direction circular variance

    dir_angle_expanded = data_angle_expanded .* (pi/180);
    ori_angle_expanded = mod(data_angle_expanded,180)*2*pi/180;

    dir_vector = abs(sum(vrsp .* exp(1i .* dir_angle_expanded),5) ./ sum(vrsp,5));
    ori_vector = abs(sum(vrsp .* exp(1i .* ori_angle_expanded),5) ./ sum(vrsp,5));
    dir_cv = 1-dir_vector;
    ori_cv = 1-ori_vector;
    dir_cv_value(:,:,:,:,rsp) = dir_cv;
    ori_cv_value(:,:,:,:,rsp) = ori_cv;

    clear vrsp noise_sigma prsp multiprsp;

end

[rp_expanded,alpha_expanded,sig_expanded,rsp_expanded] = ndgrid(gpuArray(grid_size.Rp),gpuArray(grid_size.Alpha),gpuArray(grid_size.Sig),gpuArray(grid_size.Rsp));

E = exp(-0.5*180^2./sig_expanded.^2);
E2 = exp(-0.5*90^2./sig_expanded.^2);
R_theta_pref = rsp_expanded+rp_expanded.*(1+alpha_expanded.*E);
R_theta_pref180 = rsp_expanded+rp_expanded.*(E+alpha_expanded);
R_theta_pref90 = rsp_expanded+rp_expanded.*(E2+alpha_expanded.*E2);

% direction index
% di_numeriator = Resp(theta_pref) - Resp(theta_pref+180)
di_numerator = R_theta_pref - R_theta_pref180;
% di_denominator = Resp(theta_pref)
%       = Rsp + Rp + Alpha*Rp*E
di_denominator = R_theta_pref;
di = min(di_numerator./di_denominator,1);

%orientation index
% oi_numerator = Resp(theta_pref)+Resp(theta_pref_180)-(Resp(theta_pref-90)+Resp(theta_pref-90))
oi_numerator = R_theta_pref+R_theta_pref180 - 2*R_theta_pref90;
% oi_denominator = Resp(theta_pref)+Resp(theta_pref_180)
oi_denominator = R_theta_pref+R_theta_pref180;
oi = min(oi_numerator./oi_denominator,1);

% %direction index
% di = squeeze((1 - alpha_expanded) .* (1 - E) ./ (rsp_expanded ./ rp_expanded + (1 + alpha_expanded) .* (1 + E)));
% di = min(di,1);
% 
% %orientation index
% oi = squeeze((rp_expanded.*(1+alpha_expanded).*(1+E) - 2*(rp_expanded.*(1+alpha_expanded).*E2)) ./ (4*rsp_expanded+rp_expanded.*(1 + alpha_expanded) .* (1 + E) + 2*rp_expanded.*(1+alpha_expanded).*E2));
% oi = min(oi,1);

%extract and normalize the liklihood curve of each parameter
lik_Rp = squeeze(sum(sum(sum(sum(Lik,5),4),3),2));
lik_Rp = lik_Rp./sum(lik_Rp,"all");
lik_theta_pref = squeeze(sum(sum(sum(sum(Lik,5),4),3),1))';
lik_theta_pref = lik_theta_pref./sum(lik_theta_pref,"all");
lik_Alpha = squeeze(sum(sum(sum(sum(Lik,5),4),2),1));
lik_Alpha = lik_Alpha./sum(lik_Alpha,"all");
lik_sigma = squeeze(sum(sum(sum(sum(Lik,5),3),2),1));
lik_sigma = lik_sigma./sum(lik_sigma,"all");
lik_rsp = squeeze(sum(sum(sum(sum(Lik,4),3),2),1));
lik_rsp = lik_rsp./sum(lik_rsp,"all");

assert(abs(sum(lik_Rp) - 1) < 1e-6, 'lik_Rp does not sum to 1');
assert(abs(sum(lik_theta_pref) - 1) < 1e-6, 'lik_theta_pref does not sum to 1');
assert(abs(sum(lik_Alpha) - 1) < 1e-6, 'lik_Alpha does not sum to 1');
assert(abs(sum(lik_sigma) - 1) < 1e-6, 'lik_sigma does not sum to 1');
assert(abs(sum(lik_rsp) - 1) < 1e-6, 'lik_rsp does not sum to 1');

[~,ind] = max(Lik,[],'all');
%ind = squeeze(ind);
[rp,op,alpha,sig,rsp] = ind2sub(size(Lik),ind);
ang = (0:359)';
vrsp = grid_size.Rsp(rsp) + grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)).^2./grid_size.Sig(sig).^2) + grid_size.Alpha(alpha) .* grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)+180).^2./grid_size.Sig(sig).^2);

%descriptors values in maximum likelihood condition
DI = di(rp,alpha,sig,rsp);
OI = oi(rp,alpha,sig,rsp);
DCV = dir_cv_value(rp,op,alpha,sig,rsp);
CV = ori_cv_value(rp,op,alpha,sig,rsp);


%circular variance histogram
Lik_index = sum(Lik,2);
di_bins = 0:0.05:1.0;
oi_bins = 0:0.05:1.0;
dir_cv_bins = 0:0.05:1.0;
cv_bins = 0:0.05:1.0;


di_lik = 0 * di_bins;
oi_lik = 0 * oi_bins;
dcv_lik = 0 * dir_cv_bins;
cv_lik = 0 * cv_bins;


[di_index] = discretize(di(:),di_bins);
for bin = 1:max(di_index)
    F = di_index==bin;
    di_lik(bin) = sum(Lik_index(F));
end
di_lik = di_lik./sum(di_lik);

[oi_index] = discretize(oi(:),oi_bins);
for bin = 1:max(oi_index)
    F = oi_index==bin;
    oi_lik(bin) = sum(Lik_index(F));
end
oi_lik = oi_lik./sum(oi_lik);


[dcv_index] = discretize(dir_cv_value(:),dir_cv_bins);
for bin = 1:max(dcv_index)
    F = dcv_index==bin;
    dcv_lik(bin) = sum(Lik(F));
end
dcv_lik = dcv_lik./sum(dcv_lik);

[ori_index] = discretize(ori_cv_value(:),cv_bins);
for bin = 1:max(ori_index)
    F = ori_index==bin;
    cv_lik(bin) = sum(Lik(F));
end
cv_lik = cv_lik./sum(cv_lik);

% release gpu storage
lik_theta_pref = gather(lik_theta_pref);
lik_Rp = gather(lik_Rp);
lik_rsp = gather(lik_rsp);
lik_sigma = gather(lik_sigma);
lik_Alpha = gather(lik_Alpha);
Lik = gather(Lik);
DI = gather(DI);
OI = gather(OI);
di_bins = gather(di_bins);
di_lik = gather(di_lik);
oi_bins = gather(oi_bins);
oi_lik = gather(oi_lik);
dcv_lik = gather(dcv_lik);
cv_lik = gather(cv_lik);


output_struct = struct( ...
    'noise_model',struct('type',{'proportional'},'offset',offset,'slope',slope,'k',k), ...
    'other_parameters',struct('independent_variable',{'angles'},'independent_variable_value',data.angles,...
    'mean_responses',data.mean_responses), ...
    'marginal_likelihoods',struct( ...
        'theta_pref',struct('values',grid_size.Op,'likelihoods',lik_theta_pref), ...
        'Rp',struct('values',grid_size.Rp,'likelihoods',lik_Rp), ...
        'Alpha',struct('values',grid_size.Alpha,'likelihoods',lik_Alpha), ...
        'sigma',struct('values',grid_size.Sig,'likelihoods',lik_sigma), ...
        'Rsp',struct('values',grid_size.Rsp,'likelihoods',lik_rsp)), ... 
    'maximum_likelihood_parameters',struct( ...
        'parameters',struct('theta_pref',grid_size.Op(op),'Rp',grid_size.Rp(rp),'Rn',grid_size.Alpha(alpha)*grid_size.Rp(rp),'sigma',grid_size.Sig(sig),'Rsp',grid_size.Rsp(rsp),'tunning_curve',vrsp), ...
		'descriptors',struct('di',DI,'oi',OI,'cv',CV,'dir_cv',DCV)), ...
	'descriptors',struct( ...
		'di',struct('values',di_bins,'likelihoods',di_lik), ...
		'oi',struct('values',oi_bins,'likelihoods',oi_lik), ...
		'cv',struct('values',cv_bins,'likelihoods',cv_lik), ...
		'dir_cv',struct('values',dir_cv_bins,'likelihoods',dcv_lik)));
end
