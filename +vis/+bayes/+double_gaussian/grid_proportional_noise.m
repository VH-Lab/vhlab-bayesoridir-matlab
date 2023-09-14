function [output_struct,Lik] = grid_proportional_noise(grid_size,data,noise_mdl, varargin)
% Perform a Bayesian parameter estimation for double gaussian function
% 
% BAYES_GRID_FUNCTION_PROPORTIONAL_NOISE(GRID_SIZE, DATA, NOISE_MLD)
%
% Given a set of DATA.angles (such as [0:30:360-30]) and a set of responses to each
% angle, where DATA.responses(i) is the response to a grating presented at
% DATA.angles(i), return a Bayesian estimation of model parameters to the function:
%
% R(THETA) = RSP+RP*exp(-vlt.math.angdiff(OP-ANGLES).^2/(2*SIG^2))+RN*exp(-vlt.math.angdiff(180+OP-ANGLES).^2/(2*SIG^2));
%
% Inputs:
%   GRID_SIZE: A structure with fields:
%      Rsp - the values of Rsp 
%      Rp - the values of Rp to examine
%      Alpha - the values of Rn (equal to Rp * Alpha)
%      Sig - the values of sigma
%   DATA : a structure with fields:
%      angles - the angles that were used for stimulation
%      responses - the mean responses to each stimulus
%      num_trials - the number of trials of each stimulus
%   NOISE_MLD  : Model of gaussian standard deviation: sigma = (10^(noise_mld(1)) + response^(noise_mld(2)))/sqrt(num_trials)
%   
%
% See also:
%    vis.bayes.noise.proportional()
% 

VERBOSE = 1;

vlt.data.assign(varargin{:});

% extract coefficients from noise linear regression model
offset = noise_mdl(1);
slope = noise_mdl(2);


% build bayes grid matrix

Lik = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
di = zeros(length(grid_size.Rp),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
oi = zeros(length(grid_size.Sig),1);
ori_cv_value = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
dir_cv_value = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));

% get 5 parameters' posterior probability(liklihood)
for rp = 1:length(grid_size.Rp)
    if VERBOSE,
        disp(['Loop ' int2str(rp) ' of ' int2str(length(grid_size.Rp)) '...']);
    end;
    for op = 1:length(grid_size.Op)
        for alpha = 1:length(grid_size.Alpha)
            for sig = 1:length(grid_size.Sig)
                for rsp = 1:length(grid_size.Rsp)
                    vrsp = grid_size.Rsp(rsp) + grid_size.Rp(rp) * exp(-0.5*angdiff(data.angles-grid_size.Op(op)).^2/grid_size.Sig(sig)^2) + grid_size.Alpha(alpha)* grid_size.Rp(rp) * exp(-0.5*angdiff(data.angles-(grid_size.Op(op)+180)).^2/grid_size.Sig(sig)^2);
                    prsp = normpdf(vrsp(:)-data.mean_responses,zeros(size(vrsp(:))),vis.bayes.noise.proportional(noise_mdl,vrsp(:),data.num_trials));
                    multiprsp = squeeze(prod(prsp));
                    Lik(rp,op,alpha,sig,rsp) = multiprsp;

                    %direction index
                    di(rp,alpha,sig,rsp) = (1-grid_size.Alpha(alpha)) * (1-exp(-0.5*180^2./grid_size.Sig(sig)^2))./(grid_size.Rsp(rsp)./grid_size.Rp(rp) + 1 + grid_size.Alpha(alpha)*exp(-0.5*180^2./grid_size.Sig(sig)^2));

                    %orientation index
                    oi(sig) = 1 - 2./(exp(0.5*90^2./grid_size.Sig(sig)^2)+exp(-1.5*90^2./grid_size.Sig(sig)^2));
                    %circular variance and direction circular variance
                    dir_vector = abs(sum(vrsp(:).*exp(sqrt(-1)*data.angles*pi/180))/sum(vrsp));
                    ori_vector = abs(sum(vrsp(:).*exp(sqrt(-1)*2/180*pi*mod(data.angles,180)))/sum(vrsp));
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
[M,ind] = max(Lik,[],'all');
M = squeeze(M);
ind = squeeze(ind);
[rp,op,alpha,sig,rsp] = ind2sub(size(Lik),ind);
ang = (0:359)';
vrsp = grid_size.Rsp(rsp) + grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)).^2./grid_size.Sig(sig).^2) + grid_size.Alpha(alpha) .* grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)+180).^2./grid_size.Sig(sig).^2);
%descriptors values in maximum likelihood condition
DI = di(rp,alpha,sig,rsp);
OI = oi(sig);
DCV = dir_cv_value(rp,op,alpha,sig,rsp);
CV = ori_cv_value(rp,op,alpha,sig,rsp);
di_lik = squeeze(sum(Lik,2));
%circular variance histogram
[dcv_index] = discretize(dir_cv_value(:),0:0.05:1);
for bin = 1:max(dcv_index)
    F = dcv_index==bin;
    dcv_lik(bin) = sum(Lik(F));
end
dcv_lik = dcv_lik./sum(dcv_lik);

[ori_index] = discretize(ori_cv_value(:),0:0.05:1);
for bin = 1:max(ori_index)
    F = ori_index==bin;
    ori_lik(bin) = sum(Lik(F));
end
ori_lik = ori_lik./sum(ori_lik);

output_struct = struct( ...
    'noise_model',struct('type',{'proportional'},'offset',offset,'slope',slope), ...
    'other_parameters',struct('independent_variable',{'angle'},'independent_variable_value',data.angles), ...
    'marginal_likelihood',struct( ...
        'theta_pref',struct('values',grid_size.Op,'likelihoods',lik_theta_pref), ...
        'Rp',struct('values',grid_size.Rp,'likelihoods',lik_Rp), ...
        'Rn',struct('values',grid_size.Alpha,'likelihoods',lik_Alpha), ...
        'sigma',struct('values',grid_size.Sig,'likelihoods',lik_sigma), ...
        'Rsp',struct('values',grid_size.Rsp,'likelihoods',lik_rsp)), ... 
    'maximum_likelihood',struct( ...
        'parameters',struct('theta_pref',grid_size.Op(op),'Rp',grid_size.Rp(rp),'Rn',grid_size.Alpha(alpha),'sigma',grid_size.Sig(sig),'Rsp',grid_size.Rsp(rsp),'tunning_curve',vrsp), ...
        'descriptors',struct('di',DI,'oi',OI,'cv',CV,'dir_cv',DCV)), ...
    'descriptors',struct( ...
        'di',struct('values',di,'likelihoods',di_lik), ...
        'oi',struct('values',oi,'likelihoods',lik_sigma), ...
        'cv',struct('values',ori_cv_value,'likelihoods',ori_lik), ...
        'dir_cv',struct('values',dir_cv_value,'likelihoods',dcv_lik)));
end
