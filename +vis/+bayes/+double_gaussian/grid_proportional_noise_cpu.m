function [output_struct,Lik] = grid_proportional_noise_cpu(grid_size,data,noise_mdl, options)
% Perform a Bayesian parameter estimation for double gaussian function
% 
% GRID_PROPORTIONAL_NOISE_CPU(GRID_SIZE, DATA, NOISE_MLD)
%
% Given a set of DATA.angles (such as [0:30:360-30]) and a set of responses to each
% angle, where DATA.responses(i) is the response to a grating presented at
% DATA.angles(i), return a Bayesian estimation of model parameters to the function:
%
% R(THETA) = RSP+RP*exp(-vlt.math.angdiff(OP-ANGLES).^2/(2*SIG^2))+RN*exp(-vlt.math.angdiff(180+OP-ANGLES).^2/(2*SIG^2));
%
% This function uses the CPU and uses no GPU functionality.
% 
% Inputs:
%   GRID_SIZE: A structure with fields:
%      Rsp - the values of Rsp 
%      Rp - the values of Rp to examine
%      Alpha - the values of Rn (equal to Rp * Alpha)
%      Op - the values of Op
%      Sig - the values of sigma'\
%      
%   DATA : a structure with fields:
%      angles - the angles that were used for stimulation
%      responses - the mean responses to each stimulus
%      num_trials - the number of trials of each stimulus
%   NOISE_MDL  : Model of gaussian standard deviation: sigma = (10^(noise_mld(1)) + response^(noise_mld(2)))/sqrt(num_trials)
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
% This function also takes an optional parameter:
%    'verbose': 0/1 should we print information about our progress? (default 1)
%
% See also:
%    vis.bayes.noise.proportional()
% 

arguments
        grid_size (1,1) struct
        data (1,1) struct
        noise_mdl double
        options.verbose (1,1) logical = true
end

% extract coefficients from noise linear regression model
offset = noise_mdl(1);
slope = noise_mdl(2);


% build bayes grid matrix

Lik = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig),length(grid_size.Rsp));
di = Lik;
oi = Lik;
ori_cv_value = Lik;
dir_cv_value = Lik;

% get 5 parameters' posterior probability(liklihood)
for rp = 1:length(grid_size.Rp)
    if options.verbose,
        disp(['Loop ' int2str(rp) ' of ' int2str(length(grid_size.Rp)) '...']);
    end;
    for op = 1:length(grid_size.Op)
        for alpha = 1:length(grid_size.Alpha)
            for sig = 1:length(grid_size.Sig)
                E = exp(-0.5*180^2./grid_size.Sig(sig)^2);
                E2 = exp(-0.5*90^2./grid_size.Sig(sig)^2);
                for rsp = 1:length(grid_size.Rsp)
                    vrsp = grid_size.Rsp(rsp) + grid_size.Rp(rp) * exp(-0.5*angdiff(data.angles-grid_size.Op(op)).^2/grid_size.Sig(sig)^2) + grid_size.Alpha(alpha)* grid_size.Rp(rp) * exp(-0.5*angdiff(data.angles-(grid_size.Op(op)+180)).^2/grid_size.Sig(sig)^2);
                    prsp = normpdf(vrsp(:)-data.mean_responses,zeros(size(vrsp(:))),vis.bayes.noise.proportional(noise_mdl,vrsp(:),data.num_trials));
                    multiprsp = squeeze(prod(prsp));
                    Lik(rp,op,alpha,sig,rsp) = multiprsp;

                    R_theta_pref = grid_size.Rsp(rsp)+grid_size.Rp(rp)*(1+grid_size.Alpha(alpha)*E);
                    R_theta_pref180 = grid_size.Rsp(rsp)+grid_size.Rp(rp)*(E+grid_size.Alpha(alpha));
                    R_theta_pref90 = grid_size.Rsp(rsp)+grid_size.Rp(rp)*(E2+grid_size.Alpha(alpha)*E2);

                    %direction index
                    % di_numeriator = Resp(theta_pref) - Resp(theta_pref+180)
                    di_numerator = R_theta_pref - R_theta_pref180;
                    % di_denominator = Resp(theta_pref)
                    %       = Rsp + Rp + Alpha*Rp*E
                    di_denominator = R_theta_pref + R_theta_pref180;
                    di(rp,op,alpha,sig,rsp) = min(di_numerator/di_denominator,1);

                    %orientation index
                    % oi_numerator = Resp(theta_pref)+Resp(theta_pref_180)-(Resp(theta_pref-90)+Resp(theta_pref-90))
                    oi_numerator = R_theta_pref + R_theta_pref180 - 2*R_theta_pref90;
                    oi_denominator = R_theta_pref + R_theta_pref180 + 2*R_theta_pref90;
                    oi(rp,op,alpha,sig,rsp) = min(oi_numerator/oi_denominator,1);
                    
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

Lik = Lik./sum(Lik(:));

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
[rp,op,alpha,sig,rsp] = ind2sub(size(Lik),ind);
ang = (0:359)';
vrsp = grid_size.Rsp(rsp) + grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)).^2./grid_size.Sig(sig).^2) + grid_size.Alpha(alpha) .* grid_size.Rp(rp) .* exp(-0.5*angdiff(ang-grid_size.Op(op)+180).^2./grid_size.Sig(sig).^2);
%descriptors values in maximum likelihood condition
DI = di(rp,op,alpha,sig,rsp);
OI = oi(rp,op,alpha,sig,rsp);
DCV = dir_cv_value(rp,op,alpha,sig,rsp);
CV = ori_cv_value(rp,op,alpha,sig,rsp);


%circular variance histogram
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
    di_lik(bin) = sum(Lik(F));
end
di_lik = di_lik./sum(di_lik);

[oi_index] = discretize(oi(:),oi_bins);
for bin = 1:max(oi_index)
    F = oi_index==bin;
    oi_lik(bin) = sum(Lik(F));
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


output_struct = struct( ...
	'noise_model',struct('type',{'proportional'},'offset',offset,'slope',slope), ...
	'other_parameters',struct( ...
        'independent_variable',{'angle'}, ...
        'independent_variable_value',data.angles, ...
        'mean_responses',data.mean_responses), ...
	'marginal_likelihoods',struct( ...
		'theta_pref',struct('values',grid_size.Op,'likelihoods',lik_theta_pref), ...
		'Rp',struct('values',grid_size.Rp,'likelihoods',lik_Rp), ...
		'Alpha',struct('values',grid_size.Alpha,'likelihoods',lik_Alpha), ...
		'sigma',struct('values',grid_size.Sig,'likelihoods',lik_sigma), ...
		'Rsp',struct('values',grid_size.Rsp,'likelihoods',lik_rsp)), ... 
	'maximum_likelihood_parameters',struct( ...
		'parameters',struct('theta_pref',grid_size.Op(op),'Rp',grid_size.Rp(rp),'Rn',grid_size.Alpha(alpha) .* grid_size.Rp(rp),'sigma',grid_size.Sig(sig),'Rsp',grid_size.Rsp(rsp),'tunning_curve',vrsp), ...
		'descriptors',struct('di',DI,'oi',OI,'cv',CV,'dir_cv',DCV)), ...
	'descriptors',struct( ...
		'di',struct('values',di_bins,'likelihoods',di_lik), ...
		'oi',struct('values',oi_bins,'likelihoods',oi_lik), ...
		'cv',struct('values',cv_bins,'likelihoods',cv_lik), ...
		'dir_cv',struct('values',dir_cv_bins,'likelihoods',dcv_lik)));


