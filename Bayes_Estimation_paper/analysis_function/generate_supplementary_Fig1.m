clear;close all;clc;

% Generate data for Supplementary Figure 1—an example of the joint marginal
% likelihood distribution.


parameter = data.generate_fixed_parameters(10,5,180,30,1);
ang = 0:30:359;
ideal = data.generate_simulate_data(ang,parameter);
noisy = data.generate_noise(ideal,5,1);
v = noisy.responses_stddev;
m = noisy.mean_responses;
v = v(:);
m = m(:);
noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);

% Bayes Estimate Grid Range and Size
I = struct('Rp',linspace(5,15,41), ...
    'Op',170:0.5:190, ...
    'Alpha',linspace(0.3,0.9,31), ...
    'Sig',linspace(20,45,51), ...
    'Rsp',linspace(0.1,2,39));

noisy.num_trials = 5;
[bayes,Lik] = bayes_grid_function_proportional_noise_gpu(I,noisy,noise_mdl);
rp_thetapref = squeeze(sum(Lik,[3 4 5]));
rp_alpha = squeeze(sum(Lik,[2 4 5]));
rp_sigma = squeeze(sum(Lik,[2 3 5]));
thetapref_rsp = squeeze(sum(Lik,[1 3 4]));

clear Lik;
save my_supplementary_fig1_3.mat