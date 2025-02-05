clear;clc;clf;
% generate figure 5
load my_fig5_1_kcs_numtrials.mat
%%
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%%
% Noise least-square fitting
m = [];
v = [];
for i = 1:numel(data_noisy)
    for j = 1:numel(data_noisy{i})
        m = [m;data_noisy{i}(j).mean_responses(:)];
        v = [v;data_noisy{i}(j).responses_stddev(:)];
    end
end
mdl = fitlm(log10(m),log10(v)),
noise_coefficients = mdl.Coefficients{:,1};
noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);
%%
% Bayes estimation
output = [];
for i = 1:numel(data_noisy)
    for j = 1:numel(data_noisy{i})
        data_noisy{i}(j).num_trials = 5; % number of simulate trials
        [output{end+1},~] = bayes_grid_function_proportional_noise_gpu(I,data_noisy{i}(j),noise_mdl);
    end
end
clear i j;