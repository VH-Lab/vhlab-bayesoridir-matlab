clear,clc,close all;

% generate data for figure 1 - Bayesian estimation for two types of
% simulated cells with different levels of tuning

% A,B - Two example cells: one that is well-tuned, and one that is
% poorly-tuned, and least-squares fits to both with the parameters noted in
% the graph.  [Would show that for the poorly tuned cell, the uncertainty
% of the preferred direction is not captured and in fact is somewhat
% random.]

% generate ideal simulated tuning curves
parameter_w = data.generate_fixed_parameters(10,5,90,30,1);% well tuned
parameter_p = data.generate_fixed_parameters(0.1,0,90,30,1);% poorly tuned
ang = 0:30:359;
ideal_w = data.generate_simulate_data(ang,parameter_w);
ideal_p = data.generate_simulate_data(ang,parameter_p);


% generate simulated tuning curves with +50% noise
noisy_w = data.generate_noise(ideal_w,5,1);
noisy_p = data.generate_noise(ideal_p,5,1);

% Noise sigma's linear regression in logarithm space
v_w = noisy_w.responses_stddev;% well-tuned
m_w = noisy_w.mean_responses;
v_p = noisy_p.responses_stddev;% poorly-tuned
m_p = noisy_p.mean_responses;
v = [v_w(:);v_p(:)];
m = [m_w(:);m_p(:)];
noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);
%%
% C,D - Maximum likelihood Bayesian parameter estimation for the two cells
% Bayes Estimate Grid Range and Size
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,21), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%Bayes Estimation
noisy_w.num_trials = 5;
noisy_p.num_trials = 5;
tic,
[bayes_w,Lik] = bayes_grid_function_proportional_noise_gpu(I,noisy_w,noise_mdl);
[bayes_p,~] = bayes_grid_function_proportional_noise_gpu(I,noisy_p,noise_mdl);
toc,

rp_thetapref = squeeze(sum(Lik,[3 4 5]));
rp_alpha = squeeze(sum(Lik,[2 4 5]));
rp_sigma = squeeze(sum(Lik,[2 3 5]));
thetapref_rsp = squeeze(sum(Lik,[1 3 4]));

clear Lik;