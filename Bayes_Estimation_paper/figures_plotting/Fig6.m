clear all;close;clc;
% Fig 6 - von mises noise model
load my_fig6.mat
%

% using gaussian model to fit von mises model
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
title('A','FontSize',12);%von mises model