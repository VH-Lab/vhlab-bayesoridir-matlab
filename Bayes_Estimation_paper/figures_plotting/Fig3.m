clear;close all;clc;
% Fig 3 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning(OI).
%
% Example cells - Trial 5: cell 55
%                 Trial 10: cell 25
%                 Trial 15: cell 55
load my_fig3_trial15.mat
cell_type = numel(data_ideal);
measurement_num = data_noisy{1}(1).num_trials;
experiment_num = numel(data_noisy{1});
ang = data_ideal{1}.angles;
%curves parameters
rp = linspace(0,10,5);
rn = linspace(0,5,5);
alpha = zeros(size(rn));
alpha(2:end) = 0.5;
rsp = linspace(10,0,5);
angle = 45;
sigma = 30;
true_oi = 1 - (2.*rsp + 2.*(rp + rn).*exp(-0.5*90^2/30^2))./(2.*rsp + (rp + rn).*(1 + exp(-0.5*180^2/30^2)));


%%
% ESE relative to a uniform distribution
relErr_rp = zeros(experiment_num,cell_type);
relErr_alpha = relErr_rp;
relErr_rsp = relErr_rp;
relErr_angle = relErr_rp;
relErr_sigma = relErr_rp;

% centralPercentile2 prompt
% mode and central of 25%, 50% and 75%
offset = [0.25 0.5 0.75];
percentiles = [25 50 75];
count_rp_total = zeros(numel(percentiles),cell_type);
count_alpha_total = count_rp_total;
count_rsp_total = count_rp_total;
count_angle_total = count_rp_total;
count_sigma_total = count_rp_total;
% Get all pdf in output data.

for i = 1:cell_type
    for j = 1:experiment_num
        idx = (i-1)*experiment_num + j;
        pdf_rp = output(idx).marginal_likelihood.Rp.likelihoods;
        pdf_alpha = output(idx).marginal_likelihood.Alpha.likelihoods;
        pdf_rsp = output(idx).marginal_likelihood.Rsp.likelihoods;
        pdf_angle = output(idx).marginal_likelihood.theta_pref.likelihoods;
        pdf_sigma = output(idx).marginal_likelihood.sigma.likelihoods;

        grid_rp = output(idx).marginal_likelihood.Rp.values;
        grid_alpha = output(idx).marginal_likelihood.Alpha.values;
        grid_rsp = output(idx).marginal_likelihood.Rsp.values;
        grid_angle = output(idx).marginal_likelihood.theta_pref.values;
        grid_sigma = output(idx).marginal_likelihood.sigma.values;

        % ESE relative to a uniform distribution
        relErr_rp(j,i) = vis.bayes.assess.expectedErrorRelativeToUniform(grid_rp,pdf_rp,rp(i));
        relErr_alpha(j,i) = vis.bayes.assess.expectedErrorRelativeToUniform(grid_alpha,pdf_alpha,alpha(i));
        relErr_rsp(j,i) = vis.bayes.assess.expectedErrorRelativeToUniform(grid_rsp,pdf_rsp,rsp(i));
        relErr_angle(j,i) = vis.bayes.assess.expectedErrorRelativeToUniform(grid_angle,pdf_angle,angle);
        relErr_sigma(j,i) = vis.bayes.assess.expectedErrorRelativeToUniform(grid_sigma,pdf_sigma,sigma);

        count_rp = isInCentralPercentile2(grid_rp,pdf_rp,rp(i),percentiles);
        count_alpha = isInCentralPercentile2(grid_alpha,pdf_alpha,alpha(i),percentiles);
        count_rsp = isInCentralPercentile2(grid_rsp,pdf_rsp,rsp(i),percentiles);
        count_angle = checkDirectionInOrientationPercentile(grid_angle,pdf_angle,angle,percentiles);
        count_sigma = isInCentralPercentile2(grid_sigma,pdf_sigma,sigma,percentiles);
        count_rp_total(:,i) = count_rp_total(:,i) + count_rp(:);
        count_alpha_total(:,i) = count_alpha_total(:,i) + count_alpha(:);
        count_rsp_total(:,i) = count_rsp_total(:,i) + count_rsp(:);
        count_angle_total(:,i) = count_angle_total(:,i) + count_angle(:);
        count_sigma_total(:,i) = count_sigma_total(:,i) + count_sigma(:);
    end
end
% Get fraction that 'true' value sit in the central intervals
fraction_rp = count_rp_total / experiment_num;
fraction_alpha = count_alpha_total / experiment_num;
fraction_rsp = count_rsp_total / experiment_num;
fraction_angle = count_angle_total / experiment_num;
fraction_sigma = count_sigma_total / experiment_num;
% Get percentage (%) that 'true' value sit in the central intervals
percentage_rp = fraction_rp * 100;
percentage_alpha = fraction_alpha * 100;
percentage_rsp = fraction_rsp * 100;
percentage_angle = fraction_angle * 100;
percentage_sigma = fraction_sigma * 100;

%%
% Showing one example cell for each type.
curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
cell_num = 55;
for i = 1:5
figure(i),hold on; 
plot(data_ideal{i}.angles,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',1)
errorbar(data_noisy{i}(cell_num).angles,data_noisy{i}(cell_num).mean_responses,data_noisy{i}(cell_num).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output((i-1)*100+cell_num).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',1)
xlim([-5,365]),
ylim([-.5,20]),
xticks(0:100:300),
yticks(0:4:20),
set(gca,'FontSize',8)
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
plot(center,output((i-1)*100+cell_num).descriptors.oi.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
end
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8)
xlabel('OI','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('F','FontSize',12);

figure(7),
hold on,
for i = 1:5
plot(output((i-1)*100+cell_num).marginal_likelihood.Rp.values,output((i-1)*100+cell_num).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xticks(0:5:20),
yticks(0:.2:1),
set(gca,'FontSize',8)
xlabel('R_{p}','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('G','FontSize',12);

figure(8),
hold on,
for i = 1:5
plot(output((i-1)*100+cell_num).marginal_likelihood.sigma.values,output((i-1)*100+cell_num).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
yticks(0:.2:1),
set(gca,'FontSize',8)
xlabel('\sigma','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('H','FontSize',12);

figure(9),
hold on,
for i = 1:5
plot(output((i-1)*100+cell_num).marginal_likelihood.theta_pref.values,output((i-1)*100+cell_num).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
xlim([-5,365])
ylim([-0.02,1]),
yticks(0:.2:1),
set(gca,'FontSize',8)
xlabel('\theta_{pref}','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('I','FontSize',12);

figure(10),
hold on,
for i = 1:5
plot(output((i-1)*100+cell_num).marginal_likelihood.Rsp.values,output((i-1)*100+cell_num).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xticks(-10:2:10),
yticks(0:.2:1),
set(gca,'FontSize',8)
xlabel('C','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('J','FontSize',12);

interval_tag = ['25%';'50%';'75%'];

figure(11),
hold on,
b = bar(offset,percentage_alpha,'grouped', 'FaceColor', 'flat', 'EdgeColor', 'none');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
set(gca,'FontSize',8)
xlabel('Credential Intervals of \alpha');
ylabel('Fraction of True Value in the Intervals (%)');
title('K','FontSize',12);

figure(12),
hold on,
b = bar(offset,percentage_rp,'grouped', 'FaceColor', 'flat', 'EdgeColor', 'none');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
set(gca,'FontSize',8)
xlabel('Credential Intervals of R_{p}');
ylabel('Fraction of True Value in the Intervals (%)');
title('L','FontSize',12);

figure(13),
hold on,
b = bar(offset,percentage_sigma,'grouped', 'FaceColor', 'flat', 'EdgeColor', 'none');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
set(gca,'FontSize',8)
xlabel('Credential Intervals of \sigma');
ylabel('Fraction of True Value in the Intervals (%)');
title('M','FontSize',12);

figure(14),
hold on,
b = bar(offset,percentage_angle,'grouped', 'FaceColor', 'flat', 'EdgeColor', 'none');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
set(gca,'FontSize',8)
xlabel('Credential Intervals of \theta_{pref}');
ylabel('Fraction of True Value in the Intervals (%)');
title('N','FontSize',12);

figure(15),
hold on,
b = bar(offset,percentage_rsp,'grouped', 'FaceColor', 'flat', 'EdgeColor', 'none');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
set(gca,'FontSize',8)
xlabel('Credential Intervals of C');
ylabel('Fraction of True Value in the Intervals (%)');
title('O','FontSize',12);
%%
barColors = reshape(cell2mat(cellfun(@(c) sscanf(c(2:end), '%2x%2x%2x')', curvecolor, 'UniformOutput', false)), 3, [])' / 255;

mean_relErr_alpha = mean(relErr_alpha,1);
std_relErr_alpha = std(relErr_alpha,0,1);
x = 1:length(mean_relErr_alpha);

figure(16),
hold on,
b = bar(x, mean_relErr_alpha, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = barColors;
for i = 1:length(x)
    errorbar(x(i), mean_relErr_alpha(i), std_relErr_alpha(i), ...
        'Color', barColors(i,:), ...
        'LineWidth', 1.2, ...
        'CapSize', 16, ...
        'LineStyle', 'none');
end
xlim([0.5, 5.5]);
ylim([0, 1.5]);
yticks(0:0.5:1.5);
xticks(x);
xticklabels([]);
set(gca, 'FontSize', 8, 'box', 'off');
ylabel('Relative Expected Squared Error', 'FontSize', 10);
xlabel('\alpha', 'FontSize', 10);

mean_relErr_rp = mean(relErr_rp,1);
std_relErr_rp = std(relErr_rp,0,1);
x = 1:length(mean_relErr_rp);

figure(17),
hold on,
b = bar(x, mean_relErr_rp, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = barColors;
for i = 1:length(x)
    errorbar(x(i), mean_relErr_rp(i), std_relErr_rp(i), ...
        'Color', barColors(i,:), ...
        'LineWidth', 1.2, ...
        'CapSize', 16, ...
        'LineStyle', 'none');
end
xlim([0.5, 5.5]);
ylim([0, 1.5]);
yticks(0:0.5:1.5);
xticks(x);
xticklabels([]);
set(gca, 'FontSize', 8, 'box', 'off');
ylabel('Relative Expected Squared Error', 'FontSize', 10);
xlabel('R_{p}', 'FontSize', 10);

mean_relErr_sigma = mean(relErr_sigma,1);
std_relErr_sigma = std(relErr_sigma,0,1);
x = 1:length(mean_relErr_sigma);

figure(18),
hold on,
b = bar(x, mean_relErr_sigma, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = barColors;
for i = 1:length(x)
    errorbar(x(i), mean_relErr_sigma(i), std_relErr_sigma(i), ...
        'Color', barColors(i,:), ...
        'LineWidth', 1.2, ...
        'CapSize', 16, ...
        'LineStyle', 'none');
end
xlim([0.5, 5.5]);
ylim([0, 1.5]);
yticks(0:0.5:1.5);
xticks(x);
xticklabels([]);
set(gca, 'FontSize', 8, 'box', 'off');
ylabel('Relative Expected Squared Error', 'FontSize', 10);
xlabel('\sigma', 'FontSize', 10);

mean_relErr_angle = mean(relErr_angle,1);
std_relErr_angle = std(relErr_angle,0,1);
x = 1:length(mean_relErr_angle);

figure(19),
hold on,
b = bar(x, mean_relErr_angle, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = barColors;
for i = 1:length(x)
    errorbar(x(i), mean_relErr_angle(i), std_relErr_angle(i), ...
        'Color', barColors(i,:), ...
        'LineWidth', 1.2, ...
        'CapSize', 16, ...
        'LineStyle', 'none');
end
xlim([0.5, 5.5]);
ylim([0, 1.5]);
yticks(0:0.5:1.5);
xticks(x);
xticklabels([]);
set(gca, 'FontSize', 8, 'box', 'off');
ylabel('Relative Expected Squared Error', 'FontSize', 10);
xlabel('\theta_{pref}', 'FontSize', 10);

mean_relErr_rsp = mean(relErr_rsp,1);
std_relErr_rsp = std(relErr_rsp,0,1);
x = 1:length(mean_relErr_rsp);

figure(20),
hold on,
b = bar(x, mean_relErr_rsp, 0.5, 'FaceColor', 'flat', 'EdgeColor', 'none');
b.CData = barColors;
for i = 1:length(x)
    errorbar(x(i), mean_relErr_rsp(i), std_relErr_rsp(i), ...
        'Color', barColors(i,:), ...
        'LineWidth', 1.2, ...
        'CapSize', 16, ...
        'LineStyle', 'none');
end
xlim([0.5, 5.5]);
ylim([0, 1.5]);
yticks(0:0.5:1.5);
xticks(x);
xticklabels([]);
set(gca, 'FontSize', 8, 'box', 'off');
ylabel('Relative Expected Squared Error', 'FontSize', 10);
xlabel('C', 'FontSize', 10);
%%
% file fraction data
path = 'D:\GitHub\vhlab-bayesoridir-matlab\Bayes_Estimation_paper\need to check\';
doc_name = 'Fig3_trials_15.xlsx';
writematrix(fraction_rp,[path doc_name],'Sheet','rp')
writematrix(fraction_alpha,[path doc_name],'Sheet','alpha')
writematrix(fraction_rsp,[path doc_name],'Sheet','rsp')
writematrix(fraction_angle,[path doc_name],'Sheet','angle')
writematrix(fraction_sigma,[path doc_name],'Sheet','sigma')
%%
% file relErr data
path = 'D:\GitHub\vhlab-bayesoridir-matlab\Bayes_Estimation_paper\need to check\';
doc_name = 'Fig3_trials_15_relErr.xlsx';
writematrix(relErr_rp,[path doc_name],'Sheet','rp')
writematrix(relErr_alpha,[path doc_name],'Sheet','alpha')
writematrix(relErr_rsp,[path doc_name],'Sheet','rsp')
writematrix(relErr_angle,[path doc_name],'Sheet','angle')
writematrix(relErr_sigma,[path doc_name],'Sheet','sigma')
%%
%export figures
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
F3k = figure(11);
F3l = figure(12);
F3m = figure(13);
F3n = figure(14);
F3o = figure(15);
F3p = figure(16);
F3q = figure(17);
F3r = figure(18);
F3s = figure(19);
F3t = figure(20);

path = 'D:\GitHub\vhlab-bayesoridir-matlab\Bayes_Estimation_paper\figures\noise_mdl\kcs\';
suffix = 'trial15.pdf';

exportgraphics(F3a,[path 'Figure_3a_' suffix],"ContentType","vector"),
exportgraphics(F3b,[path 'Figure_3b_' suffix],"ContentType","vector"),
exportgraphics(F3c,[path 'Figure_3c_' suffix],"ContentType","vector"),
exportgraphics(F3d,[path 'Figure_3d_' suffix],"ContentType","vector"),
exportgraphics(F3e,[path 'Figure_3e_' suffix],"ContentType","vector"),
exportgraphics(F3f,[path 'Figure_3f_' suffix],"ContentType","vector"),
exportgraphics(F3g,[path 'Figure_3g_' suffix],"ContentType","vector"),
exportgraphics(F3h,[path 'Figure_3h_' suffix],"ContentType","vector"),
exportgraphics(F3i,[path 'Figure_3i_' suffix],"ContentType","vector"),
exportgraphics(F3j,[path 'Figure_3j_' suffix],"ContentType","vector"),
exportgraphics(F3k,[path 'Figure_3k_' suffix],"ContentType","vector"),
exportgraphics(F3l,[path 'Figure_3l_' suffix],"ContentType","vector"),
exportgraphics(F3m,[path 'Figure_3m_' suffix],"ContentType","vector"),
exportgraphics(F3n,[path 'Figure_3n_' suffix],"ContentType","vector"),
exportgraphics(F3o,[path 'Figure_3o_' suffix],"ContentType","vector");
exportgraphics(F3p,[path 'Figure_3p_' suffix],"ContentType","vector"),
exportgraphics(F3q,[path 'Figure_3q_' suffix],"ContentType","vector"),
exportgraphics(F3r,[path 'Figure_3r_' suffix],"ContentType","vector"),
exportgraphics(F3s,[path 'Figure_3s_' suffix],"ContentType","vector"),
exportgraphics(F3t,[path 'Figure_3t_' suffix],"ContentType","vector");