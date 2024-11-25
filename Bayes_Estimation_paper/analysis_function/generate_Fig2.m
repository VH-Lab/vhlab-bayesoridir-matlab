clear;close all;clc;

% generate data for figure 2 - Neuronal responses and Bayesian estimation
% to stimuli at different temporal frequencies in simulated cell

% A - An imagined V1 neuron that exhibits direction selectivity
% that varies with temporal frequency. On the left is shown the response of
% the cell to a low temporal frequency.

% A - Studying a direction-selective cell where direction preference changes
% with temporal frequency. (3 cells: one that is strongly responsive and
% selective to say 45 degrees; one that is slightly less responsive and
% equally selective to 45 and 225 degrees; one that is poorly responsive,
% mostly noise). Show one curve in each of 3 colors.

parameter_l = data.generate_fixed_parameters(7,1,45,30,1);% low tf
parameter_m = data.generate_fixed_parameters(5,5,45,30,1);% median tf
parameter_h = data.generate_fixed_parameters(0,0,45,30,1);% high tf
ang = 0:22.5:359;
ideal_l = data.generate_simulate_data(ang,parameter_l);
ideal_m = data.generate_simulate_data(ang,parameter_m);
ideal_h = data.generate_simulate_data(ang,parameter_h);
% generate simulated tuning curves with +50% noise
noisy_l = data.generate_noise(ideal_l,5,1);
noisy_m = data.generate_noise(ideal_m,5,1);
noisy_h = data.generate_noise(ideal_h,5,1);
%% 
% B - A row of marginal probabilities for theta_p, Rp, and DI (use color to
% indicate cell identity).A situation where Bayesian parameter estimation
% allows insight into the responses of neurons where least-squares does
% not.

% Noise least-square fitting
m = [noisy_l.mean_responses,noisy_m.mean_responses,noisy_h.mean_responses];
v = [noisy_l.responses_stddev;noisy_m.responses_stddev;noisy_h.responses_stddev];

noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);
%%
% Bayes Estimate Grid Range and Size
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%%
% Bayes Estimation
noisy_l.num_trials = 5;
noisy_m.num_trials = 5;
noisy_h.num_trials = 5;
tic,
[bayes_l,~] = bayes_grid_function_proportional_noise_gpu(I,noisy_l,noise_mdl);
[bayes_m,~] = bayes_grid_function_proportional_noise_gpu(I,noisy_m,noise_mdl);
[bayes_h,~] = bayes_grid_function_proportional_noise_gpu(I,noisy_h,noise_mdl);
toc,
