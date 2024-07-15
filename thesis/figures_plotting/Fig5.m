clear all;close;clc;
% Fig 5 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning.
cell_type = 1;
experiment_num = 10;
ang = 0:359;
%curves parameters
rp = 7.5;
rn = 3.75;
rsp = 2.5;
angle = 45;
sigma = 30;
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
true_oi = 1 - (2.*rsp + 2.*(rp + rn).*exp(-0.5*90^2/30^2))./(2.*rsp + (rp + rn).*(1 + exp(-0.5*180^2/30^2)));
true_di = (rp - rn).*(1 - exp(-0.5*180^2./sigma^2))./(rsp + rp + rn.*exp(-0.5.*180^2./sigma^2));


curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
    file_name = ['my_fig5_',int2str(i),'_new_kcs.mat'];
    load(file_name)
    figure(i),hold on;
    plot(data_ideal{1}.angle,data_ideal{1}.responses,'Color',curvecolor{i},'LineWidth',1)
    errorbar(data_noisy{1}(2).angle,data_noisy{1}(2).mean_responses,data_noisy{1}(2).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
    plot(0:359,output{2}.maximum_likelihood.parameters.tunning_curve,'k','LineWidth',1)
    xlim([-5,365]),
    ylim([-0.5,20]),
    xticks(0:100:300),
    yticks(0:4:20),
    xlabel('Direction of Stimuli Motion'),
    ylabel('Response(Hz)'),
    legend('Ideal','+50% Noise','Bayes Estimation'),
    title(titlename{i});

    figure(6),
    if exist("edges",'var'), hold on;
    plot(center,output{2}.descriptors.oi.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
    ylim([-0.01,.3]),
    yticks(0:.1:.3),
    xlabel('OI'),
    ylabel('Probability'),
    title('F');
    end

    figure(7),
    hold on,
    plot(output{2}.marginal_likelihood.Rp.values,output{2}.marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
    ylim([-0.01,.3]),
    yticks(0:.1:.3),
    xlabel('R_{p}'),
    ylabel('Probability'),
    title('G');

    figure(8),
    hold on,
    plot(output{2}.marginal_likelihood.sigma.values,output{2}.marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
    ylim([-0.01,.3]),
    yticks(0:.1:.3),
    xlabel('\sigma'),
    ylabel('Probability'),
    title('H');

    figure(9),
    hold on,
    plot(output{2}.marginal_likelihood.theta_pref.values,output{2}.marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
    xlim([-5,365]),
    ylim([-0.01,.3]),
    yticks(0:.1:.3),
    xlabel('\theta_{p}'),
    ylabel('Probability'),
    title('I');

    figure(10),
    hold on,
    plot(output{2}.marginal_likelihood.Rsp.values,output{2}.marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
    ylim([-0.01,.3]),
    yticks(0:.1:.3),
    xlabel('Response Offset'),
    ylabel('Probability'),
    title('J');
end

%%
% export figures
F5a = figure(1);
F5b = figure(2);
F5c = figure(3);
F5d = figure(4);
F5e = figure(5);
F5f = figure(6);
F5g = figure(7);
F5h = figure(8);
F5i = figure(9);
F5j = figure(10);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs\';

exportgraphics(F5a,[path 'Figure_5a_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5b,[path 'Figure_5b_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5c,[path 'Figure_5c_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5d,[path 'Figure_5d_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5e,[path 'Figure_5e_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5f,[path 'Figure_5f_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5g,[path 'Figure_5g_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5h,[path 'Figure_5h_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5i,[path 'Figure_5i_kcs.pdf'],"ContentType","vector"),
exportgraphics(F5j,[path 'Figure_5j_kcs.pdf'],"ContentType","vector"),