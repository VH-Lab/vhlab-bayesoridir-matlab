clear;close all;clc;
load my_fig1_kcs_numtrials.mat
% Fig 1 - A comparison of least-squares fitting to Bayesian parameter
% estimation

% A,B - Two example cells: one that is well-tuned, and one that is
% poorly-tuned, and least-squares fits to both with the parameters noted in
% the graph.  [Would show that for the poorly tuned cell, the uncertainty
% of the preferred direction is not captured and in fact is somewhat
% random.]
% generate ideal simulated tuning curves
% % % parameter_w = data.generate_fixed_parameters(10,5,90,30,1);% well tuned
% % % parameter_p = data.generate_fixed_parameters(1,0,90,30,1);% poorly tuned
% % % ang = 0:10:359;
% % % ideal_w = data.generate_simulate_data(ang,parameter_w);
% % % ideal_p = data.generate_simulate_data(ang,parameter_p);


% generate simulated tuning curves with +50% noise
% % % noisy_w = data.generate_noise(ideal_w,5,1);
% % % noisy_p = data.generate_noise(ideal_p,5,1);

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
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%Bayes Estimation
noisy_w.num_trials = 5;
noisy_p.num_trials = 5;
tic,
[bayes_w,~] = bayes_grid_function_proportional_noise_gpu(I,noisy_w,noise_mdl);
[bayes_p,~] = bayes_grid_function_proportional_noise_gpu(I,noisy_p,noise_mdl);
toc,

%%
figure(1);
hold on,
plot(ideal_w.angle,ideal_w.responses,'k','LineWidth',1),% show well tuned ideal curve
errorbar(noisy_w.angle,noisy_w.mean_responses,noisy_w.responses_stderr,'b*'),%show well tuned noisy points
plot(0:359,bayes_w.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',1),
xlim([-5,365]),
ylim([0,20]),
yticks(0:4:20),
set(gca,'FontSize',8),
legend('Ideal','+50% Noise','Bayes Estimation','FontSize',7),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
title('A','FontSize',12);%Well Tuned Curve

figure(2);
hold on,
plot(ideal_p.angle,ideal_p.responses,'k','LineWidth',1),% show poorly tuned ideal curve
errorbar(noisy_p.angle,noisy_p.mean_responses,noisy_p.responses_stderr,'b*'),%show poorly tuned noisy points
plot(0:359,bayes_p.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',1),
xlim([-5,365]),
ylim([0,20]),
yticks(0:4:20),
set(gca,'FontSize',8),
legend('Ideal','+50% Noise','Bayes Estimation','FontSize',7),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
title('B','FontSize',12);%Poorly Tuned Curve
%%
% E,F - Marginal likelihood of theta-pref
figure(3);
plot(bayes_w.marginal_likelihood.theta_pref.values,bayes_w.marginal_likelihood.theta_pref.likelihoods,'k','LineWidth',1)
xlim([-5,365]),
ylim([-0.02,1]),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('\theta_{pref}','FontSize',10),
ylabel('Marginal Likelihood','FontSize',10),
title('C','FontSize',12);%Well Tuned Theta-pref

figure(4);
plot(bayes_p.marginal_likelihood.theta_pref.values,bayes_p.marginal_likelihood.theta_pref.likelihoods,'k','LineWidth',1)
xlim([-5,365]),
ylim([-0.02,1]),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('\theta_{pref}','FontSize',10),
ylabel('Marginal Likelihood','FontSize',10),
title('D','FontSize',12);%Poorly Tuned Theta-pref

% G,H - Marginal likelihood of Rp
figure(5);
plot(bayes_w.marginal_likelihood.Rp.values,bayes_w.marginal_likelihood.Rp.likelihoods,'k','LineWidth',1)
ylim([-0.02,1]),
xticks(0:5:20),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('R_{p}','FontSize',10),
ylabel('Marginal Likelihood','FontSize',10),
title('E','FontSize',12);%Well Tuned Rp

figure(6);
plot(bayes_p.marginal_likelihood.Rp.values,bayes_p.marginal_likelihood.Rp.likelihoods,'k','LineWidth',1)
ylim([-0.02,1]),
xticks(0:5:20),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('R_{p}','FontSize',10),
ylabel('Marginal Likelihood','FontSize',10),
title('F','FontSize',12);%Poorly Tuned Rp

% I,J - Marginal likelihood of OI

edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;

figure(7);
plot(center,bayes_w.descriptors.oi.histogram_likelihoods,'k','LineWidth',1)
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('Orientation Index (OI)','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('G','FontSize',12);%Well Tuned OI

figure(8);
plot(center,bayes_p.descriptors.oi.histogram_likelihoods,'k','LineWidth',1)
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('Orientation Index (OI)','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('H','FontSize',12);%Poorly Tuned OI

% K,L - Marginal likelihood of DI

figure(9);
plot(center,bayes_w.descriptors.di.histogram_likelihoods,'k','LineWidth',1)
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('Direction Index (DI)','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('I','FontSize',12);%Well Tuned DI

figure(10);
plot(center,bayes_p.descriptors.di.histogram_likelihoods,'k','LineWidth',1)
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('Direction Index (DI)','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('J','FontSize',12);%Poorly Tuned DI

%%
% export figures
F1a = figure(1);
F1b = figure(2);
F1c = figure(3);
F1d = figure(4);
F1e = figure(5);
F1f = figure(6);
F1g = figure(7);
F1h = figure(8);
F1i = figure(9);
F1j = figure(10);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs\';

exportgraphics(F1a,[path 'Figure_1a_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1b,[path 'Figure_1b_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1c,[path 'Figure_1c_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1d,[path 'Figure_1d_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1e,[path 'Figure_1e_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1f,[path 'Figure_1f_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1g,[path 'Figure_1g_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1h,[path 'Figure_1h_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1i,[path 'Figure_1i_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F1j,[path 'Figure_1j_kcs_numtrials.pdf'],"ContentType","vector")
