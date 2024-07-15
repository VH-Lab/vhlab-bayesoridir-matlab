% Reset workspace and load exist data documents from the file
clear all;close;clc;
%load my_fig7_debug_3_before.mat;
load my_fig7_new_kcs.mat
load stevesolddata.mat;
if exist("noisy_BME","var")
    noisy_before = noisy_BME;
end
%%
% Fig 7 - Bayes Estimation Result from Real Neuron Dataset
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

good_both = intersect(good_indexes_before,good_indexes_motionexposure);
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
ylim([-0.1,0.3]),
yticks(-0.1:0.1:0.3)
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
yticks(0:.2:1),
xlabel('OI'),
ylabel('Probability'),
title('F');

figure(7),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.Rp.values,output(i).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([-.01 .31])
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('R_{p}'),
ylabel('Probability'),
title('G');

figure(8),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.sigma.values,output(i).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([0 61])
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('\sigma'),
ylabel('Probability'),
title('H');

figure(9),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.theta_pref.values,output(i).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([-5 365])
ylim([-0.02,1]),
yticks(0:.2:1),
xlabel('\theta_{p}'),
ylabel('Probability'),
title('I');

figure(10),
hold on,
for i = 1:5
plot(output(i).marginal_likelihood.Rsp.values,output(i).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([-.11 .21])
ylim([-0.02,1]),
yticks(0:.2:1),
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
for i=1:numel(good_both),
    for j=1:numel(cell_structures{good_both(i)}),
        if strcmp(cell_structures{good_both(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params'),
            data_bootsbefore{i} = cell_structures{good_both(i)}(j).data;
        end;
    end;
end;
for i = 1:numel(data_bootsbefore),
    rsp_min(i) = min(data_bootsbefore{i}(:,1));
    rsp_max(i) = max(data_bootsbefore{i}(:,1));
    rp_min(i) = min(data_bootsbefore{i}(:,2));
    rp_max(i) = max(data_bootsbefore{i}(:,2));
    op_min(i) = min(data_bootsbefore{i}(:,3));
    op_max(i) = max(data_bootsbefore{i}(:,3));
    sig_min(i) = min(data_bootsbefore{i}(:,4));
    sig_max(i) = max(data_bootsbefore{i}(:,4));
    rn_min(i) = min(data_bootsbefore{i}(:,5));
    rn_max(i) = max(data_bootsbefore{i}(:,5));
end;
rsp_range = [min(rsp_min),max(rsp_max)]
rp_range = [min(rp_min),max(rp_max)]
op_range = [min(op_min),max(op_max)]
sig_range = [min(sig_min),max(sig_max)]
rn_range = [min(rn_min),max(rn_max)]

var = max(noisy_before(1).mean_responses);

N = histc(data_bootsbefore{3}(:,1),linspace(-var,var,40));
N = N./sum(N,'all');
figure(11),
bar(linspace(-var,var,40),N);
ylim([0,1]),
xlabel('Response Offset'),
ylabel('Probability'),

N = histc(data_bootsbefore{3}(:,2),linspace(0.001,3*var,60));
N = N./sum(N,'all');
figure(12),
bar(linspace(0.001,3*var,60),N);
ylim([0,1]),
xlabel('R_{p}'),
ylabel('Probability'),

N = histc(data_bootsbefore{3}(:,3),0:5:359);
N = N./sum(N,'all');
figure(13),
bar(0:5:359,N);
ylim([0,1]),
xlabel('\theta_{p}'),
ylabel('Probability'),

N = histc(data_bootsbefore{3}(:,4),linspace(1,60,60));
N = N./sum(N,'all');
figure(14),
bar(linspace(1,60,60),N);
ylim([0,1]),
xlabel('\sigma'),
ylabel('Probability'),

N = histc(data_bootsbefore{3}(:,5),linspace(0,3*var,60));
N = N./sum(N,'all');
figure(15),
bar(linspace(0,3*var,60),N);
ylim([0,1]),
xlabel('R_{null}'),
ylabel('Probability'),

%%
% export figures
F7a = figure(1);
F7b = figure(2);
F7c = figure(3);
F7d = figure(4);
F7e = figure(5);
F7f = figure(6);
F7g = figure(7);
F7h = figure(8);
F7i = figure(9);
F7j = figure(10);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs\';

exportgraphics(F7a,[path 'Figure_7a_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7b,[path 'Figure_7b_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7c,[path 'Figure_7c_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7d,[path 'Figure_7d_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7e,[path 'Figure_7e_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7f,[path 'Figure_7f_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7g,[path 'Figure_7g_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7h,[path 'Figure_7h_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7i,[path 'Figure_7i_kcs.pdf'],"ContentType","vector"),
exportgraphics(F7j,[path 'Figure_7j_kcs.pdf'],"ContentType","vector"),