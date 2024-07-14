clear all;close;clc;
load my_fig4_new.mat;
% Fig 4 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning(DI).

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
plot(data_ideal{i}.angle,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',1)
errorbar(data_mean{i}(1).angle,data_mean{i}(1).mean_responses,data_mean{i}(1).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*10-9).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',1)
xlim([-5,365]),
ylim([-.5,20]),
xticks(0:100:300),
yticks(0:4:20),
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
plot(center,output(i*10-3).descriptors.di.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
end
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('DI'),
ylabel('Probability'),
title('F');

figure(7),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rp.values,output(i*10-3).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('R_{p}'),
ylabel('Probability'),
title('G');

figure(8),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.sigma.values,output(i*10-3).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('\sigma'),
ylabel('Probability'),
title('H');

figure(9),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.theta_pref.values,output(i*10-3).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([-5,365]),
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('\theta_{p}'),
ylabel('Probability'),
title('I');

figure(10),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rsp.values,output(i*10-3).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('Response Offset'),
ylabel('Probability'),
title('J');
%%
% export figures
F4a = figure(1);
F4b = figure(2);
F4c = figure(3);
F4d = figure(4);
F4e = figure(5);
F4f = figure(6);
F4g = figure(7);
F4h = figure(8);
F4i = figure(9);
F4j = figure(10);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs';

exportgraphics(F4a,[path 'Figure_4a.pdf'],"ContentType","vector"),
exportgraphics(F4b,[path 'Figure_4b.pdf'],"ContentType","vector"),
exportgraphics(F4c,[path 'Figure_4c.pdf'],"ContentType","vector"),
exportgraphics(F4d,[path 'Figure_4d.pdf'],"ContentType","vector"),
exportgraphics(F4e,[path 'Figure_4e.pdf'],"ContentType","vector"),
exportgraphics(F4f,[path 'Figure_4f.pdf'],"ContentType","vector"),
exportgraphics(F4g,[path 'Figure_4g.pdf'],"ContentType","vector"),
exportgraphics(F4h,[path 'Figure_4h.pdf'],"ContentType","vector"),
exportgraphics(F4i,[path 'Figure_4i.pdf'],"ContentType","vector"),
exportgraphics(F4j,[path 'Figure_4j.pdf'],"ContentType","vector"),