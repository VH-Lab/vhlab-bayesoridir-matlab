clear;close;clc;
% Fig 3 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning(OI).
cell_type = 5;
measurement_num = 5;
experiment_num = 10;
ang = 0:45:359;
rp = linspace(0,10,5);
rn = linspace(0,5,5);
rsp = linspace(10,0,5);
for i = 1:5
parameters1 = data.generate_fixed_parameters(rp(i),rn(i),45,30,rsp(i));
end
true_oi = 1 - (2.*rsp + 2.*(rp + rn).*exp(-0.5*90^2/30^2))./(2.*rsp + (rp + rn).*(1 + exp(-0.5*180^2/30^2)));

load my_fig3_new.mat

curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
plot(data_ideal{i}.angle,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',1)
errorbar(data_mean{i}(7).angle,data_mean{i}(7).mean_responses,data_mean{i}(7).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*10-3).maximum_likelihood.parameters.tunning_curve,'k')
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
plot(center,output(i*10-3).descriptors.oi.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
end
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('OI'),
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
xlim([-5,365])
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
F3a = figure(1);
F3b = figure(2);
F3c = figure(3);
F3d = figure(4);
F3e = figure(5);
F3f = figure(6);
F3g = figure(7);
F3h = figure(8);
F3i = figure(9);
F3j = figure(10);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\';

exportgraphics(F3a,[path 'Figure_3a.pdf'],"ContentType","vector"),
exportgraphics(F3b,[path 'Figure_3b.pdf'],"ContentType","vector"),
exportgraphics(F3c,[path 'Figure_3c.pdf'],"ContentType","vector"),
exportgraphics(F3d,[path 'Figure_3d.pdf'],"ContentType","vector"),
exportgraphics(F3e,[path 'Figure_3e.pdf'],"ContentType","vector"),
exportgraphics(F3f,[path 'Figure_3f.pdf'],"ContentType","vector"),
exportgraphics(F3g,[path 'Figure_3g.pdf'],"ContentType","vector"),
exportgraphics(F3h,[path 'Figure_3h.pdf'],"ContentType","vector"),
exportgraphics(F3i,[path 'Figure_3i.pdf'],"ContentType","vector"),
exportgraphics(F3j,[path 'Figure_3j.pdf'],"ContentType","vector"),