% Reset workspace and load exist data documents from the file
clear all;close;clc;
load my_fig7_debug_2.mat;
load stevesolddata.mat;
%%
% Fig 4 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning.
good_indexes_before = [];
good_indexes_motionexposure = [];

for i=1:numel(cell_structures),
    for j=1:numel(cell_structures{i}),
        if strcmp(cell_structures{i}(j).type,'TP Ach OT vec varies p'),
            if cell_structures{i}(j).data<0.05,
                good_indexes_before(end+1) = i;
            end;
        end;
        if strcmp(cell_structures{i}(j).type,'TP ME Ach OT vec varies p'),
            if cell_structures{i}(j).data<0.05,
                good_indexes_motionexposure(end+1) = i;
            end;
        end;
    end;
end;
% create plotting edges(n+1) and center(n) N = curves number
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;

% extract di likelihood data in output structure. Transform the data form into
% plotting suitable form

curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
errorbar(noisy_before(i).angle,noisy_before(i).mean_responses,noisy_before(i).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i).maximum_likelihood.parameters.tunning_curve,'Color',curvecolor{i},'LineWidth',1)
xlim([-5,365]),
ylim([-0.12,0.2]),
yticks(0:0.05:0.2)
xlabel('Direction of Stimuli Motion'),
ylabel('Response(Hz)'),
legend('V1 Data Points','Bayes Estimation'),
title(titlename{i});
end
%%
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
figure(6),hold on;
for i = 1:5
plot(center,output(i).descriptors.oi.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
end
ylim([-0.02,1]),
xlabel('OI'),
ylabel('Probability'),
title('F');

figure(7),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.Rp.values,output(i).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xlabel('R_{p}'),
ylabel('Probability'),
title('G');

figure(8),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.sigma.values,output(i).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xlabel('\sigma'),
ylabel('Probability'),
title('H');

figure(9),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.theta_pref.values,output(i).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xlabel('\theta_{p}'),
ylabel('Probability'),
title('I');

figure(10),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.Rsp.values,output(i).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xlabel('Response Offset'),
ylabel('Probability'),
title('J');

%%
% To see a list of all response data available

all_types = {cell_structures{1}(:).type}';


% Example here

% data column 1: Response Offset
% column 2 : Rpref
% column 3 : theta_pref
% column 4 : sigma (tuning width)
% column 5 : Rnull
data_bootsbefore = [];
for i=1:numel(good_indexes_before),
    for j=1:numel(cell_structures{good_indexes_before(i)}),
        if strcmp(cell_structures{good_indexes_before(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params'),
            data_bootsbefore{i} = cell_structures{good_indexes_before(i)}(j).data;
        end;
    end;
end;

N = histc(data_bootsbefore{3}(:,1),linspace(0.001,0.3,30));
% N = N./sum(N,'all');
figure(11),
bar(linspace(0.001,0.3,30),N);
