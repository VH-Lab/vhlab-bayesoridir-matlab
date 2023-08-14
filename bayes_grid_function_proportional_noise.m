function [output_struct,Lik] = bayes_grid_function_proportional_noise(grid_size,data,noise_mdl)
% extract coefficients of noise linear regression model
offset = noise_mdl(1);
slope = noise_mdl(2);
% build bayes grid matrix
Lik = zeros(length(grid_size.Rp),length(grid_size.Op),length(grid_size.Alpha),length(grid_size.Sig));

% get 4 parameters' posterior probability(liklihood)
for i = 1:length(grid_size.Rp)
    i
    for j = 1:length(grid_size.Op)
        for k = 1:length(grid_size.Alpha)
            for h = 1:length(grid_size.Sig)
                vrsp = grid_size.Rp(i) * exp(-0.5*angdiff(data.angle-grid_size.Op(j)).^2/grid_size.Sig(h)^2) + grid_size.Alpha(k)* grid_size.Rp(i) * exp(-0.5*angdiff(data.angle-(grid_size.Op(j)+180)).^2/grid_size.Sig(h)^2);
                prsp = normpdf(vrsp-data.mean_responses',zeros(size(vrsp)),10.^offset*vrsp.^slope);
                multiprsp = squeeze(prod(prsp));
                Lik(i,j,k,h) = multiprsp;
            end
        end
    end
end

%extract and normalize the liklihood curve of each parameter
lik_Rp = squeeze(sum(sum(sum(Lik,4),3),2));
lik_Rp = lik_Rp./sum(lik_Rp,"all");
lik_theta_pref = squeeze(sum(sum(sum(Lik,4),3),1))';
lik_theta_pref = lik_theta_pref./sum(lik_theta_pref,"all");
lik_Alpha = squeeze(sum(sum(sum(Lik,4),2),1));
lik_Alpha = lik_Alpha./sum(lik_Alpha,"all");
lik_sigma = squeeze(sum(sum(sum(Lik,3),2),1));
lik_sigma = lik_sigma./sum(lik_sigma,"all");

[M,ind] = max(Lik,[],'all');
M = squeeze(M);
ind = squeeze(ind);
[i,j,k,h] = ind2sub(size(Lik),ind);
% size(i),size(j),size(k),size(h),
ang = (0:359)';
vrsp = grid_size.Rp(i) .* exp(-0.5*angdiff(ang-grid_size.Op(j)).^2./grid_size.Sig(h).^2) + grid_size.Alpha(k) .* grid_size.Rp(i) .* exp(-0.5*angdiff(ang-grid_size.Op(j)+180).^2./grid_size.Sig(h).^2);
% size(lik.vrsp)
DI = (1 - grid_size.Alpha);
LDI = lik_Alpha;

output_struct = struct( ...
    'noise_model',struct('type',{'proportional'},'offset',offset,'slope',slope), ...
    'other_parameters',struct('independent_variable',{'angle'},'independent_variable_value',data.angle), ...
    'marginal_likelihood',struct( ...
        'theta_pref',struct('values',grid_size.Op,'likelihoods',lik_theta_pref), ...
        'Rp',struct('values',grid_size.Rp,'likelihoods',lik_Rp), ...
        'Rn',struct('values',grid_size.Alpha,'likelihoods',lik_Alpha), ...
        'sigma',struct('values',grid_size.Sig,'likelihoods',lik_sigma), ...
        'Rsp',struct('values',data.mean_responses,'likelihoods',vrsp)), ... %response_mean have 8 values and vrsp have 360
    'maximum_likelihood_parameters',struct('theta_pref',grid_size.Op(j),'Rp',grid_size.Rp(i),'Rn',grid_size.Alpha(k),'sigma',grid_size.Sig(h),'Rsp',max(vrsp)), ...
    'descriptors',struct( ...
        'di',struct('values',DI,'likelihoods',LDI), ...
        'oi',struct('values',[],'likelihoods',[]), ...
        'cv',struct('values',[],'likelihoods',[]), ...
        'dir_cv',struct('values',[],'likelihoods',[])));
end