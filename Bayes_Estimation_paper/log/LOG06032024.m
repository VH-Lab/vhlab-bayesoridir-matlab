clear,clc,clf;
% DATE（MMDDYYYY） - 06052024
% Create GPU core function - bayes_grid_function_proportional_noise_gpu
% Goal - When fitting grid size is too large and there are huge amount of
% samples, GPU acceleration can compress time into a acceptable range
% instead of a month. Compare two outputs between previous and gpu verson.

% Solution - use gpuArray and change the expression of fitting_curve
% equation in the loop.
%          - Check Lik distribution to see if there is any difference.

%extract experiment data
load stevesolddata.mat
good_indexes_before = [];
good_indexes_motionexposure = [];

for i=1:numel(cell_structures),
    for j=1:numel(cell_structures{i}),
        if strcmp(cell_structures{i}(j).type,'TP Ach OT vec varies p'),
            if cell_structures{i}(j).data<0.05,
                good_indexes_before(end+1) = i;
            end;
        end;
        if strcmp(cell_structures{i}(j).type,'TP ME Ach OT vec varies p'),
            if cell_structures{i}(j).data<0.05,
                good_indexes_motionexposure(end+1) = i;
            end;
        end;
    end;
end;

good_both = intersect(good_indexes_before,good_indexes_motionexposure);

% To extract direction tuning data for each cell before motion exposure, use:
data_before = [];
for i=1:numel(good_both),
    for j=1:numel(cell_structures{good_both(i)}),
        if strcmp(cell_structures{good_both(i)}(j).type,'TP Ach OT Response curve'),
            data_before{end+1} = cell_structures{good_both(i)}(j).data;
            noisy_before(i).angle = data_before{end}(1,:);
            noisy_before(i).mean_responses = data_before{end}(2,:)';
            noisy_before(i).responses_stddev = data_before{end}(3,:);
            noisy_before(i).responses_stderr = data_before{end}(4,:);
        end;
    end;
end;

% Bayesian Estimation
% bayes grid input
var = max(noisy_before(3).mean_responses);

I = struct('Rp',linspace(0.001,3*var,30), ...
    'Op',0:10:359, ...
    'Alpha',linspace(0,1,21), ...
    'Sig',linspace(1,60,30), ...
    'Rsp',linspace(-var,var,20));

%noise fitting model
m = [];
v = [];
for i = 1:numel(noisy_before)
    v = [v,noisy_before(i).responses_stderr];
    m = [m,noisy_before(i).mean_responses];
end
v = v(m>0);
m = m(m>0);
mdl = fitlm(log10(m),log10(v)),
figure(),plot(mdl);
xlabel('log10(response mean)'),
ylabel('log10(response stddev)')
noise_coefficients = mdl.Coefficients{:,1};
%%
[output2,lik2] = bayes_grid_function_proportional_noise_gpu(I,noisy_before(3),noise_coefficients);
[output1,lik1] = bayes_grid_function_proportional_noise(I,noisy_before(3),noise_coefficients);

size(lik1),
size(lik2)
figure
plot(lik1(1:10),lik2(1:10))
%%
% likelihood comparison
figure(1),hold on
plot(output1.marginal_likelihood.Alpha.values,output1.marginal_likelihood.Alpha.likelihoods)
plot(output2.marginal_likelihood.Alpha.values,output2.marginal_likelihood.Alpha.likelihoods+0.01)

figure(2),hold on
plot(output1.marginal_likelihood.Rp.values,output1.marginal_likelihood.Rp.likelihoods)
plot(output2.marginal_likelihood.Rp.values,output2.marginal_likelihood.Rp.likelihoods+0.01)

figure(3),hold on
plot(output1.marginal_likelihood.sigma.values,output1.marginal_likelihood.sigma.likelihoods)
plot(output2.marginal_likelihood.sigma.values,output2.marginal_likelihood.sigma.likelihoods+0.01)

figure(4),hold on
plot(output1.marginal_likelihood.theta_pref.values,output1.marginal_likelihood.theta_pref.likelihoods)
plot(output2.marginal_likelihood.theta_pref.values,output2.marginal_likelihood.theta_pref.likelihoods+0.01)

figure(5),hold on
plot(output1.marginal_likelihood.Rsp.values,output1.marginal_likelihood.Rsp.likelihoods)
plot(output2.marginal_likelihood.Rsp.values,output2.marginal_likelihood.Rsp.likelihoods+0.01)
