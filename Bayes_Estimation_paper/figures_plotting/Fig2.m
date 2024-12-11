clear;close all;clc;
load my_fig2_kcs_numtrials.mat
% Fig 2 - Bayesian parameter estimation allows the rigorous study of tuning
% curve properties across different conditions, even if the cells exhibit
% poor tuning

% A - An imagined V1 neuron that exhibits direction selectivity
% that varies with temporal frequency. On the left is shown the response of
% the cell to a low temporal frequency.

% A - Studying a direction-selective cell where direction preference changes
% with temporal frequency. (3 cells: one that is strongly responsive and
% selective to say 45 degrees; one that is slightly less responsive and
% equally selective to 45 and 225 degrees; one that is poorly responsive,
% mostly noise). Show one curve in each of 3 colors.

figure(1);
hold on,
plot(0:359,bayes_l.maximum_likelihood.parameters.tunning_curve,'Color','#0072BD','LineWidth',1),% show low tf bayes curve
errorbar(noisy_l.angle,noisy_l.mean_responses,noisy_l.responses_stderr,'*','Color','#A2142F','MarkerSize',7),
xlim([-5,365]),
ylim([0,10]),
yticks(0:2:10),
set(gca,'FontSize',8),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
legend('Bayes','+50% Noise','FontSize',7),
title('A','FontSize',12);% Low TF Curve
hold off;

figure(2);
hold on,
plot(0:359,bayes_m.maximum_likelihood.parameters.tunning_curve,'Color','#EDB120','LineWidth',1)% show median tf bayes curve
errorbar(noisy_m.angle,noisy_m.mean_responses,noisy_m.responses_stderr,'*','Color','#A2142F','MarkerSize',7)
xlim([-5,365]),
ylim([0,10]),
yticks(0:2:10),
set(gca,'FontSize',8),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
legend('Bayes','+50% Noise','FontSize',7)
title('B','FontSize',12);% Median TF Curve
hold off;

figure(3);
hold on,
plot(0:359,bayes_h.maximum_likelihood.parameters.tunning_curve,'Color','#77AC30','LineWidth',1),% show high tf bayes curve
errorbar(noisy_h.angle,noisy_h.mean_responses,noisy_h.responses_stddev,'*','Color','#A2142F','MarkerSize',7),
xlim([-5,365])
ylim([0,10]),
yticks(0:2:10),
set(gca,'FontSize',8),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
legend('Bayes','+50% Noise','FontSize',7),
title('C','FontSize',12);% High TF Curve
hold off;


%% 
figure(4);
hold on,
plot(bayes_l.marginal_likelihood.theta_pref.values,bayes_l.marginal_likelihood.theta_pref.likelihoods,'Color','#0072BD','LineWidth',1),
plot(bayes_m.marginal_likelihood.theta_pref.values,bayes_m.marginal_likelihood.theta_pref.likelihoods,'Color','#EDB120','LineWidth',1),
plot(bayes_h.marginal_likelihood.theta_pref.values,bayes_h.marginal_likelihood.theta_pref.likelihoods,'Color','#77AC30','LineWidth',1),
xlim([-5,365])
ylim([-0.02,.7]),
yticks(0:.1:.7),
set(gca,'FontSize',8),
xlabel('\theta_{pref}','FontSize',10),
ylabel('Marginal Likelihood','FontSize',10),
legend('Low TF','Median TF','High TF','FontSize',7),
title('D','FontSize',12); % Marginal Probability of \theta_{p}
hold off;

figure(5);
hold on,
plot(bayes_l.marginal_likelihood.Rp.values,bayes_l.marginal_likelihood.Rp.likelihoods,'Color','#0072BD','LineWidth',1),
plot(bayes_m.marginal_likelihood.Rp.values,bayes_m.marginal_likelihood.Rp.likelihoods,'Color','#EDB120','LineWidth',1),
plot(bayes_h.marginal_likelihood.Rp.values,bayes_h.marginal_likelihood.Rp.likelihoods,'Color','#77AC30','LineWidth',1),
ylim([-0.02,.5]),
xticks(0:5:20),
yticks(0:.1:.5),
set(gca,'FontSize',8),
xlabel('R_{p}','FontSize',10),
ylabel('Marginal Likelihood','FontSize',10),
legend('Low TF','Median TF','High TF','FontSize',7),
title('E','FontSize',12); % Marginal Probability of R_{p}
hold off;

edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;

figure(6);
hold on,
plot(center,bayes_l.descriptors.di.histogram_likelihoods,'Color','#0072BD','LineWidth',1),
plot(center,bayes_m.descriptors.di.histogram_likelihoods,'Color','#EDB120','LineWidth',1),
plot(center,bayes_h.descriptors.di.histogram_likelihoods,'Color','#77AC30','LineWidth',1),
% plot(center,bayes_l_di,'Color','#0072BD','LineWidth',2);
% plot(center,bayes_m_di,'Color','#EDB120','LineWidth',2);
% plot(center,bayes_h_di,'Color','#77AC30','LineWidth',2);
ylim([-0.02,.5]),
xticks(0:.2:1),
yticks(0:.1:.5),
set(gca,'FontSize',8),
xlabel('DI','FontSize',10),
ylabel('Likelihood','FontSize',10),
legend('Low TF','Median TF','High TF','FontSize',7),
title('F','FontSize',12); % Marginal Probability of DI
hold off;

%% 
% export figures
F2a = figure(1);
F2b = figure(2);
F2c = figure(3);
F2d = figure(4);
F2e = figure(5);
F2f = figure(6);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs\';

exportgraphics(F2a,[path 'Figure_2a_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F2b,[path 'Figure_2b_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F2c,[path 'Figure_2c_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F2d,[path 'Figure_2d_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F2e,[path 'Figure_2e_kcs_numtrials.pdf'],"ContentType","vector")
exportgraphics(F2f,[path 'Figure_2f_kcs_numtrials.pdf'],"ContentType","vector")