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
Lik = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
di = zeros(length(grid_size.Rp),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
oi = zeros(length(grid_size.Rp),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
ori_cv_value = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
dir_cv_value = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
% get 5 parameters' posterior probability(liklihood)
for rp = 1:length(grid_size.Rp)
    fprintf('This loop round left: %d \n',length(grid_size.Rp)+1-rp);
    for op = 1:length(grid_size.Op)
        for alpha = 1:length(grid_size.Alpha)
            for sig = 1:length(grid_size.Sig)
                for rsp = 1:length(grid_size.Rsp)
                    fitting_rsp_v = grid_size.Rsp(rsp) + grid_size.Rp(rp) * exp(-0.5*angdiff(data.angle-grid_size.Op(op)).^2/grid_size.Sig(sig)^2) + grid_size.Alpha(alpha)* grid_size.Rp(rp) * exp(-0.5*angdiff(data.angle-(grid_size.Op(op)+180)).^2/grid_size.Sig(sig)^2);
                    prsp = normpdf(data.mean_responses',fitting_rsp_v,10.^offset*fitting_rsp_v.^slope);
                    multiprsp = squeeze(prod(prsp));
                    Lik(rp,op,alpha,sig,rsp) = multiprsp;
                    %direction index
                    di(rp,alpha,sig,rsp) = (1-grid_size.Alpha(alpha)) * (1-exp(-0.5*180^2./grid_size.Sig(sig)^2))./(grid_size.Rsp(rsp)./grid_size.Rp(rp) + 1 + grid_size.Alpha(alpha)*exp(-0.5*180^2./grid_size.Sig(sig)^2));
                    %orientation index
                    oi(rp,alpha,sig,rsp) = 1 - (2*grid_size.Rsp(rsp) + 2*grid_size.Rp(rp)*(1+grid_size.Alpha(alpha))*exp(-0.5*90^2/grid_size.Sig(sig)^2))./(2*grid_size.Rsp(rsp) + grid_size.Rp(rp)*(1+grid_size.Alpha(alpha))*(1 + exp(-0.5*90^2/grid_size.Sig(sig)^2)));
                    %circular variance and direction circular variance
                    dir_angle = (0:359)';
                    ori_angle = [0:179,0:179]';
                    dir_vrsp = grid_size.Rsp(rsp) + grid_size.Rp(rp) * exp(-0.5*angdiff(dir_angle-grid_size.Op(op)).^2/grid_size.Sig(sig)^2) + grid_size.Alpha(alpha)* grid_size.Rp(rp) * exp(-0.5*angdiff(dir_angle-(grid_size.Op(op)+180)).^2/grid_size.Sig(sig)^2);
                    dir_angle = dir_angle.*(pi/180);
                    dir_vector = abs(sum(dir_vrsp.*exp(1i.*dir_angle))/sum(dir_vrsp));
                    ori_vector = abs(sum(dir_vrsp.*exp(1i.*ori_angle))/sum(dir_vrsp));
                    dir_cv = 1-dir_vector;
                    ori_cv = 1-ori_vector;
                    dir_cv_value(rp,op,alpha,sig,rsp) = dir_cv;
                    ori_cv_value(rp,op,alpha,sig,rsp) = ori_cv;
                end
            end
        end
    end
end
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
[~,ind] = max(Lik,[],'all');
ind = squeeze(ind);
[rp,op,alpha,sig,rsp] = ind2sub(size(Lik),ind);
ang = (0:359)';
fitting_rsp_v = grid_size.Rsp(rsp) + grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)).^2./grid_size.Sig(sig).^2) + grid_size.Alpha(alpha) .* grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)+180).^2./grid_size.Sig(sig).^2);

%descriptors values in maximum likelihood condition
DI = di(rp,alpha,sig,rsp);
OI = oi(rp,alpha,sig,rsp);
DCV = dir_cv_value(rp,op,alpha,sig,rsp);
CV = ori_cv_value(rp,op,alpha,sig,rsp);
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

[dcv_index] = discretize(dir_cv_value(:),0:0.05:1);
for bin = 1:max(dcv_index)
    F = dcv_index==bin;
    dcv_count(bin) = sum(F);
    dcv_lik(bin) = sum(Lik(F));
end
dcv_lik = dcv_lik./sum(dcv_lik);

[cv_index] = discretize(ori_cv_value(:),0:0.05:1);
for bin = 1:max(cv_index)
    F = cv_index==bin;
    cv_count(bin) = sum(F);
    cv_lik(bin) = sum(Lik(F));
end
cv_lik = cv_lik./sum(cv_lik);

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
        'descriptors',struct('di',DI,'oi',OI,'cv',CV,'dir_cv',DCV)), ...
    'descriptors',struct( ...
        'di',struct('histogram_count',di_count,'histogram_likelihoods',di_lik_count), ...
        'oi',struct('histogram_count',oi_count,'histogram_likelihoods',oi_lik_count), ...
        'cv',struct('histogram_count',cv_count,'histogram_likelihoods',cv_lik), ...
        'dir_cv',struct('histogram_count',dcv_count,'histogram_likelihoods',dcv_lik)));
end
