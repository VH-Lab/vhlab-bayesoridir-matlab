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
I = struct('Rp',linspace(0.1,20,30), ...
    'Op',0:10:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(5,60,30), ...
    'Rsp',linspace(0.1,10,20));
%Bayes Estimation
tic,
[bayes_w,lik_w] = bayes_grid_function_proportional_noise(I,noisy_w,noise_coefficients);
[bayes_p,lik_p] = bayes_grid_function_proportional_noise(I,noisy_p,noise_coefficients);
toc,
%%
% C,D - Maximum Likelihood Curve
clear all;close;clc;
load Figure1.mat
F1a = figure(1);
hold on,
plot(ideal_w.angle,ideal_w.responses,'k','LineWidth',2),% show well tuned ideal curve
errorbar(noisy_w.angle,noisy_w.mean_responses,noisy_w.responses_stderr,'b*'),%show well tuned noisy points
plot(0:359,bayes_w.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',2),
xlim([0,370]),
ylim([0,20]),
legend('Ideal','+50% Noise','Bayes Estimation'),
xlabel('Direction of Stimuli Motion'),
ylabel('Response(Hz)'),
title('A');%Well Tuned Curve
f1a = gca;

F1b = figure(2);
hold on,
plot(ideal_p.angle,ideal_p.responses,'k','LineWidth',2),% show poorly tuned ideal curve
errorbar(noisy_p.angle,noisy_p.mean_responses,noisy_p.responses_stderr,'b*'),%show poorly tuned noisy points
plot(0:359,bayes_p.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',2),
xlim([0,370]),
ylim([0,20]),
legend('Ideal','+50% Noise','Bayes Estimation'),
xlabel('Direction of Stimuli Motion'),
ylabel('Response(Hz)'),
title('B');%Poorly Tuned Curve
f1b = gca;
%%
% E,F - Marginal likelihood of theta-pref
F1c = figure(3);
plot(bayes_w.marginal_likelihood.theta_pref.values,bayes_w.marginal_likelihood.theta_pref.likelihoods,'LineWidth',2)
ylim([0,1]),
xlabel('\theta_{pref}'),
ylabel('Marginal Likelihood'),
title('C');%Well Tuned Theta-pref
f1c = gca;

F1d = figure(4);
plot(bayes_p.marginal_likelihood.theta_pref.values,bayes_p.marginal_likelihood.theta_pref.likelihoods,'LineWidth',2)
ylim([0,1])
xlabel('\theta_{pref}'),
ylabel('Marginal Likelihood'),
title('D');%Poorly Tuned Theta-pref
f1d = gca;

% G,H - Marginal likelihood of Rp
F1e = figure(5);
plot(bayes_w.marginal_likelihood.Rp.values,bayes_w.marginal_likelihood.Rp.likelihoods,'LineWidth',2)
ylim([0,1])
xlabel('R_{pref}'),
ylabel('Marginal Likelihood'),
title('E');%Well Tuned Rp
f1e = gca;

F1f = figure(6);
plot(bayes_p.marginal_likelihood.Rp.values,bayes_p.marginal_likelihood.Rp.likelihoods,'LineWidth',2)
ylim([0,1])
xlabel('R_{pref}'),
ylabel('Marginal Likelihood'),
title('F');%Poorly Tuned Rp
f1f = gca;

% I,J - Marginal likelihood of OI
edges = 0:0.01:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
F1g = figure(7);
bar(center,bayes_w.descriptors.oi.histogram_likelihoods,1)
ylim([0,0.5])
xlabel('Orientation Index (OI)'),
ylabel('Probability'),
title('G');%Well Tuned OI
f1g = gca;

F1h = figure(8);
bar(center,bayes_p.descriptors.oi.histogram_likelihoods,1)
ylim([0,0.5])
xlabel('Orientation Index (OI)'),
ylabel('Probability'),
title('H');%Poorly Tuned OI
f1h = gca;

% K,L - Marginal likelihood of DI
F1i = figure(9);
bar(center,bayes_w.descriptors.di.histogram_likelihoods,1)
ylim([0,0.5])
xlabel('Direction Index (DI)'),
ylabel('Probability'),
title('I');%Well Tuned DI
f1i = gca;

F1j = figure(10);
bar(center,bayes_p.descriptors.di.histogram_likelihoods,1)
ylim([0,0.5])
xlabel('Direction Index (DI)'),
ylabel('Probability'),
title('J');%Poorly Tuned DI
f1j = gca;

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
%% 
% All in 1 figure
fnew = figure(200);
F1a_copy = copyobj(f1a,fnew);
F1b_copy = copyobj(f1b,fnew);
F1c_copy = copyobj(f1c,fnew);
F1d_copy = copyobj(f1d,fnew);
F1e_copy = copyobj(f1e,fnew);
F1f_copy = copyobj(f1f,fnew);
F1g_copy = copyobj(f1g,fnew);
F1h_copy = copyobj(f1h,fnew);
F1i_copy = copyobj(f1i,fnew);
F1j_copy = copyobj(f1j,fnew);

subplot(5,2,1,F1a_copy);
legend('Ideal','+50% Noise','Bayes Estimation'),
subplot(5,2,2,F1b_copy);
legend('Ideal','+50% Noise','Bayes Estimation'),
subplot(5,2,3,F1c_copy);
subplot(5,2,4,F1d_copy);
subplot(5,2,5,F1e_copy);
subplot(5,2,6,F1f_copy);
subplot(5,2,7,F1g_copy);
subplot(5,2,8,F1h_copy);
subplot(5,2,9,F1i_copy);
subplot(5,2,10,F1j_copy);
