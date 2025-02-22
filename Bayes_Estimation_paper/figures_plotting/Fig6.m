clear all;close;clc;
% Fig 6 - von mises noise model
load my_fig6.mat
%
cell_num = 5;
experiment_num = 10;
ang = 0:45:359;
%curves parameters
rp = linspace(10,0,cell_num);
angle = 45;
sigma = 1;
for i = 1:cell_num
    parameters = data.vonmises.generate_von_mises_fixed_parameters(rp(i),angle,sigma);
    data_ideal{i} = data.vonmises.generate_von_mises_data(0:359,parameters);
end
%plotting
curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
plot(data_ideal{i}.angle,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',1)
errorbar(data_noisy{i}(6).angle,data_noisy{i}(6).mean_responses,data_noisy{i}(6).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*10-4).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',1)
xlim([-5,365]),
ylim([-.5,20]),
xticks(0:100:300),
yticks(0:4:20),
set(gca,'FontSize',8),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
legend('Ideal','+50% Noise','Bayes Estimation','FontSize',7),
title(titlename{i},'FontSize',12);
end
%%
% create plotting edges(n+1) and center(n) N = curves number
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;

figure(6),hold on;
for i = 1:5
plot(center,output(i*10-4).descriptors.oi.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
end
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('OI','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('F','FontSize',12);

figure(7),
hold on,
for i = 1:5
plot(output(i*10-4).marginal_likelihood.Rp.values,output(i*10-4).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xticks(0:2.5:10), 
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('R_{p}','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('G','FontSize',12);

figure(8),
hold on,
for i = 1:5
plot(output(i*10-4).marginal_likelihood.sigma.values,output(i*10-4).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('\sigma','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('H','FontSize',12);

figure(9),
hold on,
for i = 1:5
plot(output(i*10-4).marginal_likelihood.theta_pref.values,output(i*10-4).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([-5,365]),
ylim([-0.02,1]),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('\theta_{pref}','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('I','FontSize',12);

figure(10),
hold on,
for i = 1:5
plot(output(i*10-4).marginal_likelihood.Rsp.values,output(i*10-4).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xticks(0:2:10),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('C','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('J','FontSize',12);
%%
% export figures
% F6a = figure(1);
% F6b = figure(2);
% F6c = figure(3);
% F6d = figure(4);
% F6e = figure(5);
% F6f = figure(6);
% F6g = figure(7);
% F6h = figure(8);
% F6i = figure(9);
% F6j = figure(10);
% 
% path = 'D:\GitHub\vhlab-bayesoridir-matlab\Bayes_Estimation_paper\figures\noise_mdl\kcs\';
% 
% exportgraphics(F6a,[path 'Figure_6a_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6b,[path 'Figure_6b_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6c,[path 'Figure_6c_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6d,[path 'Figure_6d_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6e,[path 'Figure_6e_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6f,[path 'Figure_6f_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6g,[path 'Figure_6g_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6h,[path 'Figure_6h_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6i,[path 'Figure_6i_kcs_numtrials.pdf'],"ContentType","vector")
% exportgraphics(F6j,[path 'Figure_6j_kcs_numtrials.pdf'],"ContentType","vector")