clear;close;clc;

% Supplementary Figure 1 — An example of the joint marginal likelihood
% distribution.
load my_supplementary_fig1.mat

figure(1);
hold on,
plot(ideal.angle,ideal.responses,'k','LineWidth',1),% show well tuned ideal curve
errorbar(noisy.angle,noisy.mean_responses,noisy.responses_stderr,'b*'),%show well tuned noisy points
plot(0:359,bayes.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',1),
xlim([-5,365]),
ylim([0,20]),
yticks(0:4:20),
set(gca,'FontSize',8),
legend('Ideal','+50% Noise','Bayes Estimation','FontSize',7),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
title('A','FontSize',12);

figure(2);
h = heatmap(bayes.marginal_likelihood.theta_pref.values,bayes.marginal_likelihood.Rp.values,rp_thetapref);
h.XLabel = 'Marginal Likelihood of \theta_{pref}';
h.YLabel = 'Marginal Likelihood of R_{p}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';

figure(3);
h = heatmap(bayes.marginal_likelihood.Alpha.values,bayes.marginal_likelihood.Rp.values,rp_alpha);
h.XLabel = 'Marginal Likelihood of Alpha';
h.YLabel = 'Marginal Likelihood of R_{p}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';

figure(4);
h = heatmap(bayes.marginal_likelihood.sigma.values,bayes.marginal_likelihood.Rp.values,rp_sigma);
h.XLabel = 'Marginal Likelihood Sigma';
h.YLabel = 'Marginal Likelihood R_{p}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';

figure(5);
h = heatmap(bayes.marginal_likelihood.Rsp.values,bayes.marginal_likelihood.theta_pref.values,thetapref_rsp);
h.XLabel = 'Marginal Likelihood of Rsp';
h.YLabel = 'Marginal Likelihood of \theta_{pref}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';
