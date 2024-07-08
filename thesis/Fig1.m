clear all;close;clc;
% Fig 1 - A comparison of least-squares fitting to Bayesian parameter
% estimation

% A,B - Two example cells: one that is well-tuned, and one that is
% poorly-tuned, and least-squares fits to both with the parameters noted in
% the graph.  [Would show that for the poorly tuned cell, the uncertainty
% of the preferred direction is not captured and in fact is somewhat
% random.]
% generate ideal simulated tuning curves
parameter_w = data.generate_fixed_parameters(10,5,90,30,1);% well tuned
parameter_p = data.generate_fixed_parameters(1,0,90,30,1);% poorly tuned
ang = 0:10:359;
ideal_w = data.generate_simulate_data(ang,parameter_w);
ideal_p = data.generate_simulate_data(ang,parameter_p);



% generate simulated tuning curves with +50% noise
noisy_w = data.generate_noise(ideal_w,5,1);
noisy_p = data.generate_noise(ideal_p,5,1);

% Noise sigma's linear regression in logarithm space
v_w = [];% well-tuned
m_w = [];
v_w = noisy_w.responses_stddev;
m_w = noisy_w.mean_responses;
v_p = [];% poorly-tuned
m_p = [];
v_p = noisy_p.responses_stddev;
m_p = noisy_p.mean_responses;
v = [v_w;v_p];
m = [m_w m_p];
mdl = fitlm(log10(m),log10(v)),
figure(100),plot(mdl);
xlabel('log10(response mean)'),
ylabel('log10(response stddev)'),
title('Poorly Tuned'),
noise_coefficients = mdl.Coefficients{:,1};
%%
% C,D - Maximum likelihood Bayesian parameter estimation for the two cells
% Bayes Estimate Grid Range and Size
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%Bayes Estimation
tic,
[bayes_w,lik_w] = bayes_grid_function_proportional_noise_gpu(I,noisy_w,noise_coefficients);
[bayes_p,lik_p] = bayes_grid_function_proportional_noise_gpu(I,noisy_p,noise_coefficients);
toc,

%%
% C,D - Maximum Likelihood Curve
clear all;close;clc;
load my_fig1.mat

%%
F1a = figure(1);
hold on,
plot(ideal_w.angle,ideal_w.responses,'k','LineWidth',1),% show well tuned ideal curve
errorbar(noisy_w.angle,noisy_w.mean_responses,noisy_w.responses_stderr,'b*'),%show well tuned noisy points
plot(0:359,bayes_w.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',1),
xlim([-5,365]),
ylim([0,20]),
legend('Ideal','+50% Noise','Bayes Estimation'),
xlabel('Direction of Stimuli Motion'),
ylabel('Response(Hz)'),
title('A');%Well Tuned Curve

F1b = figure(2);
hold on,
plot(ideal_p.angle,ideal_p.responses,'k','LineWidth',1),% show poorly tuned ideal curve
errorbar(noisy_p.angle,noisy_p.mean_responses,noisy_p.responses_stderr,'b*'),%show poorly tuned noisy points
plot(0:359,bayes_p.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',1),
xlim([-5,365]),
ylim([0,20]),
legend('Ideal','+50% Noise','Bayes Estimation'),
xlabel('Direction of Stimuli Motion'),
ylabel('Response(Hz)'),
title('B');%Poorly Tuned Curve
%%
% E,F - Marginal likelihood of theta-pref
F1c = figure(3);
plot(bayes_w.marginal_likelihood.theta_pref.values,bayes_w.marginal_likelihood.theta_pref.likelihoods,'k')
xlim([-5,365]),
ylim([-0.02,1]),
xlabel('\theta_{pref}'),
ylabel('Marginal Likelihood'),
title('C');%Well Tuned Theta-pref

F1d = figure(4);
plot(bayes_p.marginal_likelihood.theta_pref.values,bayes_p.marginal_likelihood.theta_pref.likelihoods,'k')
xlim([-5,365]),
ylim([-0.02,1])
xlabel('\theta_{pref}'),
ylabel('Marginal Likelihood'),
title('D');%Poorly Tuned Theta-pref

% G,H - Marginal likelihood of Rp
F1e = figure(5);
plot(bayes_w.marginal_likelihood.Rp.values,bayes_w.marginal_likelihood.Rp.likelihoods,'k')
ylim([-0.02,1])
xlabel('R_{pref}'),
ylabel('Marginal Likelihood'),
title('E');%Well Tuned Rp

F1f = figure(6);
plot(bayes_p.marginal_likelihood.Rp.values,bayes_p.marginal_likelihood.Rp.likelihoods,'k')
ylim([-0.02,1])
xlabel('R_{pref}'),
ylabel('Marginal Likelihood'),
title('F');%Poorly Tuned Rp

% I,J - Marginal likelihood of OI

edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
F1g = figure(7);
plot(center,bayes_w.descriptors.oi.histogram_likelihoods,'k')
ylim([-0.02,1])
xlabel('Orientation Index (OI)'),
ylabel('Probability'),
title('G');%Well Tuned OI

F1h = figure(8);
plot(center,bayes_p.descriptors.oi.histogram_likelihoods,'k')
ylim([-0.02,1])
xlabel('Orientation Index (OI)'),
ylabel('Probability'),
title('H');%Poorly Tuned OI

% K,L - Marginal likelihood of DI

F1i = figure(9);
plot(center,bayes_w.descriptors.di.histogram_likelihoods,'k')
ylim([-0.02,1])
xlabel('Direction Index (DI)'),
ylabel('Probability'),
title('I');%Well Tuned DI

F1j = figure(10);
plot(center,bayes_p.descriptors.di.histogram_likelihoods,'k')
ylim([-0.02,1])
xlabel('Direction Index (DI)'),
ylabel('Probability'),
title('J');%Poorly Tuned DI

%%
% export figures
path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\';

exportgraphics(F1a,[path 'Figure_1a.pdf'],"ContentType","vector")
exportgraphics(F1b,[path 'Figure_1b.pdf'],"ContentType","vector")
exportgraphics(F1c,[path 'Figure_1c.pdf'],"ContentType","vector")
exportgraphics(F1d,[path 'Figure_1d.pdf'],"ContentType","vector")
exportgraphics(F1e,[path 'Figure_1e.pdf'],"ContentType","vector")
exportgraphics(F1f,[path 'Figure_1f.pdf'],"ContentType","vector")
exportgraphics(F1g,[path 'Figure_1g.pdf'],"ContentType","vector")
exportgraphics(F1h,[path 'Figure_1h.pdf'],"ContentType","vector")
exportgraphics(F1i,[path 'Figure_1i.pdf'],"ContentType","vector")
exportgraphics(F1j,[path 'Figure_1j.pdf'],"ContentType","vector")
