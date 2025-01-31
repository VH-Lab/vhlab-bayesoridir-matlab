clear,clc,close all;



parameter = data.generate_von_mises_fixed_parameters(10,pi,0.5);
ang = 0:30:359;
ideal = data.generate_von_mises_data(ang,parameter);

noisy = data.generate_von_mises_noise(ideal,5,1);

v = noisy.responses_stddev;
m = noisy.mean_responses;
noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);

% Bayes Estimate Grid Range and Size
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,30), ...
    'Rsp',linspace(0.1,10,30));

noisy.num_trials = 5;

[bayes,~] = bayes_grid_function_proportional_noise_gpu(I,noisy,noise_mdl);