function [output_struct,Lik] = bayes_grid_function_proportional_noise(grid_size,data,noise_mdl)
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
offset = noise_mdl(1);
slope = noise_mdl(2);
% build bayes grid matrix
Lik = gpuArray(zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp)));
di = gpuArray(zeros(length(grid_size.Rp),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp)));
oi = gpuArray(zeros(length(grid_size.Rp),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp)));
% ori_cv_value = gpuArray(zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp)));
% dir_cv_value = gpuArray(zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp)));
% get 5 parameters' posterior probability(liklihood)
[rp_expanded,op_expanded,alpha_expanded,sig_expanded] = ndgrid(gpuArray(grid_size.Rp),gpuArray(grid_size.Op),gpuArray(grid_size.Alpha),gpuArray(grid_size.Sig));
data_angle_expanded = gpuArray(reshape(data.angle, 1, 1, 1, 1, []));
data_mean_expanded = gpuArray(reshape(data.mean_responses', 1, 1, 1, 1, []));
for rsp = 1:length(grid_size.Rsp)
    fprintf('This loop round left: %d \n',length(grid_size.Rsp)+1-rsp);
    fitting_rsp_v = grid_size.Rsp(rsp) + rp_expanded .* exp(-0.5 * angdiff(data_angle_expanded - op_expanded) .^ 2 ./ sig_expanded .^ 2) + alpha_expanded .* rp_expanded .* exp(-0.5 * angdiff(data_angle_expanded - (op_expanded + 180)) .^ 2 ./ sig_expanded .^ 2);
    prsp = normpdf(data_mean_expanded,fitting_rsp_v,10.^offset*abs(fitting_rsp_v).^slope);
    multiprsp = squeeze(prod(prsp,5));
    Lik(:,:,:,:,rsp) = multiprsp;

    %circular variance and direction circular variance
    % dir_angle = (0:5:359)';
    % dir_angle_expanded = gpuArray(reshape(dir_angle, 1, 1, 1, 1, 1, []));
    % ori_angle = [0:5:179,0:5:179]';
    % ori_angle_expanded = gpuArray(reshape(ori_angle, 1, 1, 1, 1, 1, []));
    % dir_vrsp = grid_size.Rsp(rsp) + rp_expanded .* exp(-0.5 * angdiff(dir_angle_expanded - op_expanded) .^ 2 ./ sig_expanded .^ 2) + alpha_expanded .* rp_expanded .* exp(-0.5 * angdiff(dir_angle_expanded - (op_expanded + 180)) .^ 2 ./ sig_expanded .^ 2);
    % dir_angle_expanded = dir_angle_expanded .* (pi/180);
    % dir_vector = abs(sum(dir_vrsp .* exp(1i .* dir_angle_expanded)) ./ sum(dir_vrsp(:)));
    % ori_vector = abs(sum(dir_vrsp .* exp(1i .* ori_angle_expanded)) ./ sum(dir_vrsp(:)));
    % dir_cv = 1-dir_vector;
    % ori_cv = 1-ori_vector;
    % dir_cv_value(:,:,:,:,rsp) = dir_cv;
    % ori_cv_value(:,:,:,:,rsp) = ori_cv;
end

[rp_expanded,alpha_expanded,sig_expanded,rsp_expanded] = ndgrid(gpuArray(grid_size.Rp),gpuArray(grid_size.Alpha),gpuArray(grid_size.Sig),gpuArray(grid_size.Rsp));
%direction index
di = squeeze((1 - alpha_expanded) .* (1 - exp(-0.5 * 180 ^ 2 ./ sig_expanded .^ 2)) ./ (rsp_expanded ./ rp_expanded + 1 + alpha_expanded .* exp(-0.5 * 180 ^ 2 ./ sig_expanded .^ 2)));
%orientation index
oi = squeeze(1 - (2 .* rsp_expanded + 2 .* rp_expanded .* (1 + alpha_expanded) .* exp(-0.5 * 90 ^ 2 ./ sig_expanded .^ 2)) ./ (2 * rsp_expanded + rp_expanded .* (1 + alpha_expanded) .* (1 + exp(-0.5 * 180 ^ 2 ./ sig_expanded .^ 2))));

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
ind = squeeze(ind);
[rp,op,alpha,sig,rsp] = ind2sub(size(Lik),ind);
ang = (0:359)';
fitting_rsp_v = grid_size.Rsp(rsp) + grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)).^2./grid_size.Sig(sig).^2) + grid_size.Alpha(alpha) .* grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)+180).^2./grid_size.Sig(sig).^2);

%descriptors values in maximum likelihood condition
DI = di(rp,alpha,sig,rsp);
OI = oi(rp,alpha,sig,rsp);

% DCV = dir_cv_value(rp,op,alpha,sig,rsp);
% CV = ori_cv_value(rp,op,alpha,sig,rsp);
di_lik = squeeze(sum(Lik,2));
oi_lik = squeeze(sum(Lik,2));

%circular variance histogram
[di_index] = discretize(di(:),0:0.05:1);
for bin = 1:max(di_index)
    F = di_index==bin;
    di_count(bin) = sum(F);
    di_lik_count(bin) = sum(di_lik(F));
end
di_lik_count = di_lik_count./sum(di_lik_count);

[oi_index] = discretize(oi(:),0:0.05:1);
for bin = 1:max(oi_index)
    F = oi_index==bin;
    oi_count(bin) = sum(F);
    oi_lik_count(bin) = sum(oi_lik(F));
end
oi_lik_count = oi_lik_count./sum(oi_lik_count);

% [dcv_index] = discretize(dir_cv_value(:),0:0.05:1);
% for bin = 1:max(dcv_index)
%     F = dcv_index==bin;
%     dcv_count(bin) = sum(F);
%     dcv_lik(bin) = sum(Lik(F));
% end
% dcv_lik = dcv_lik./sum(dcv_lik);
% 
% [cv_index] = discretize(ori_cv_value(:),0:0.05:1);
% for bin = 1:max(cv_index)
%     F = cv_index==bin;
%     cv_count(bin) = sum(F);
%     cv_lik(bin) = sum(Lik(F));
% end
% cv_lik = cv_lik./sum(cv_lik);

% release gpu storage
lik_theta_pref = gather(lik_theta_pref);
lik_Rp = gather(lik_Rp);
lik_rsp = gather(lik_rsp);
lik_sigma = gather(lik_sigma);
lik_Alpha = gather(lik_Alpha);
Lik = gather(Lik);
DI = gather(DI);
OI = gather(OI);
di_count = gather(di_count);
di_lik_count = gather(di_lik_count);
oi_count = gather(oi_count);
oi_lik_count = gather(oi_lik_count);

output_struct = struct( ...
    'noise_model',struct('type',{'proportional'},'offset',offset,'slope',slope), ...
    'other_parameters',struct('independent_variable',{'angle'},'independent_variable_value',data.angle), ...
    'marginal_likelihood',struct( ...
        'theta_pref',struct('values',grid_size.Op,'likelihoods',lik_theta_pref), ...
        'Rp',struct('values',grid_size.Rp,'likelihoods',lik_Rp), ...
        'Alpha',struct('values',grid_size.Alpha,'likelihoods',lik_Alpha), ...
        'sigma',struct('values',grid_size.Sig,'likelihoods',lik_sigma), ...
        'Rsp',struct('values',grid_size.Rsp,'likelihoods',lik_rsp)), ... 
    'maximum_likelihood',struct( ...
        'parameters',struct('theta_pref',grid_size.Op(op),'Rp',grid_size.Rp(rp),'Rn',grid_size.Alpha(alpha)*grid_size.Rp(rp),'sigma',grid_size.Sig(sig),'Rsp',grid_size.Rsp(rsp),'tunning_curve',fitting_rsp_v), ...
        'descriptors',struct('di',DI,'oi',OI)), ...
    'descriptors',struct( ...
        'di',struct('histogram_count',di_count,'histogram_likelihoods',di_lik_count), ...
        'oi',struct('histogram_count',oi_count,'histogram_likelihoods',oi_lik_count)));
end
