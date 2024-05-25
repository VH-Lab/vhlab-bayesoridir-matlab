clear all;close;clc;
% Fig 3 -

ang = 0:45:359;
rp = linspace(0,10,5);
rn = linspace(0,5,5);
rsp = linspace(10,0,5);
for i = 1:5
parameters1 = data.generate_fixed_parameters(rp(i),rn(i),45,30,rsp(i));
end
true_oi = 1 - (2.*rsp + 2.*(rp + rn).*exp(-0.5*90^2/30^2))./(2.*rsp + (rp + rn).*(1 + exp(-0.5*180^2/30^2)));

load my_fig3.mat

curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
plot(data_ideal{i}.angle,data_ideal{i}.responses,'Color',curvecolor{i})
errorbar(data_mean{i}(7).angle,data_mean{i}(7).mean_responses,data_mean{i}(7).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*10-3).maximum_likelihood.parameters.tunning_curve,'k')
xlim([0,370]),
ylim([0,20]),
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
bayes_oi{i} = sum(reshape(output(i*10-3).descriptors.oi.histogram_likelihoods,5,[]));
plot(center,bayes_oi{i},'Color',curvecolor{i},'LineWidth',2)
end
ylim([0,1]),
xlabel('OI'),
ylabel('Probability'),
title('F');

figure(7),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rp.values,output(i*10-3).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([0,1]),
xlabel('R_{p}'),
ylabel('Probability'),
title('G');

figure(8),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.sigma.values,output(i*10-3).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([0,1]),
xlabel('\sigma'),
ylabel('Probability'),
title('H');

figure(9),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.theta_pref.values,output(i*10-3).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([0,1]),
xlabel('\theta_{p}'),
ylabel('Probability'),
title('I');

figure(10),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rsp.values,output(i*10-3).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',2),
end
ylim([0,1]),
xlabel('Response Offset'),
ylabel('Probability'),
title('J');