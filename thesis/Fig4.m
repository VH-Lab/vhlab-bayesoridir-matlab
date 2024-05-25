clear all;close;clc;
load my_fig4.mat;
% Fig 4 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning.

cell_num = 5;
experiment_num = 10;
ang = 0:45:359;
%curves parameters
rp = 10;
rn = linspace(10,0,cell_num);
rsp = 0;
angle = 45;
sigma = 30;
% create plotting edges(n+1) and center(n) N = curves number
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
% calculate the 'true DI value' from all 50 simulated curves.
true_di = (rp - rn).*(1 - exp(-0.5*180^2./sigma^2))./(rsp + rp + rn.*exp(-0.5.*180^2./sigma^2));
% extract di likelihood data in output structure. Transform the data form into
% plotting suitable form
for i = 1:5
parameters1 = data.generate_fixed_parameters(rp,rn(i),angle,sigma,rsp);
end

curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
plot(data_ideal{i}.angle,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',2)
errorbar(data_mean{i}(1).angle,data_mean{i}(1).mean_responses,data_mean{i}(1).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*10-9).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',2)
xlim([-5,365]),
ylim([-0.5,20]),
xlabel('Direction of Stimuli Motion'),
ylabel('Response(Hz)'),
legend('Ideal','+50% Noise','Bayes Estimation'),
title(titlename{i});
end
%%
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
figure(6),hold on;
for i = 1:5
bayes_di{i} = sum(reshape(output(i*10-3).descriptors.di.histogram_likelihoods,5,[]));
plot(center,bayes_di{i},'Color',curvecolor{i},'LineWidth',2)
end
ylim([-0.02,1]),
xlabel('DI'),
ylabel('Probability'),
title('F');

figure(7),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rp.values,output(i*10-3).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([-0.02,1]),
xlabel('R_{p}'),
ylabel('Probability'),
title('G');

figure(8),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.sigma.values,output(i*10-3).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([-0.02,1]),
xlabel('\sigma'),
ylabel('Probability'),
title('H');

figure(9),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.theta_pref.values,output(i*10-3).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([-0.02,1]),
xlabel('\theta_{p}'),
ylabel('Probability'),
title('I');

figure(10),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rsp.values,output(i*10-3).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([-0.02,1]),
xlabel('Response Offset'),
ylabel('Probability'),
title('J');