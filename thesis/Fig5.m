clear all;close;clc;
% Fig 4 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning.
cell_type = 1;
experiment_num = 1;
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
    file_name = ['Figure5_',int2str(i),'.mat'];
    load(file_name)
    figure(i),hold on;
    plot(data_ideal{1}.angle,data_ideal{1}.responses,'Color',curvecolor{i},'LineWidth',2)
    errorbar(data_noisy{1}(2).angle,data_noisy{1}(2).mean_responses,data_noisy{1}(2).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
    plot(0:359,output(2).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',2)
    xlim([-5,365]),
    ylim([-0.5,20]),
    xlabel('Direction of Stimuli Motion'),
    ylabel('Response(Hz)'),
    legend('Ideal','+50% Noise','Bayes Estimation'),
    title(titlename{i});

    figure(6),
    if exist("edges"), hold on;
    bayes_oi = sum(reshape(output(2).descriptors.di.histogram_likelihoods,5,[]));
    plot(center,bayes_oi,'Color',curvecolor{i},'LineWidth',2)
    ylim([-0.02,1]),
    xlabel('OI'),
    ylabel('Probability'),
    title('F');
    end
    figure(7),
    hold on,
    plot(output(2).marginal_likelihood.Rp.values,output(2).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',2),
    ylim([-0.02,1]),
    xlabel('R_{p}'),
    ylabel('Probability'),
    title('G');

    figure(8),
    hold on,
    plot(output(2).marginal_likelihood.sigma.values,output(2).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',2),
    ylim([-0.02,1]),
    xlabel('\sigma'),
    ylabel('Probability'),
    title('H');

    figure(9),
    hold on,
    plot(output(2).marginal_likelihood.theta_pref.values,output(2).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',2),
    ylim([-0.02,1]),
    xlabel('\theta_{p}'),
    ylabel('Probability'),
    title('I');

    figure(10),
    hold on,
    plot(output(2).marginal_likelihood.Rsp.values,output(2).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',2),
    ylim([-0.02,1]),
    xlabel('Response Offset'),
    ylabel('Probability'),
    title('J');
end

