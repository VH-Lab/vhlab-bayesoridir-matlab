clear;close all;clc;
% Fig 3 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning(OI).
cell_type = 5;
measurement_num = 5;
experiment_num = 100;
ang = 0:45:359;
%curves parameters
rp = linspace(0,10,5);
rn = linspace(0,5,5);
alpha = zeros(size(rn));
alpha(2:end) = 0.5;
rsp = linspace(10,0,5);
angle = 45;
sigma = 30;
for i = 1:5
parameters1 = data.generate_fixed_parameters(rp(i),rn(i),45,30,rsp(i));
end
true_oi = 1 - (2.*rsp + 2.*(rp + rn).*exp(-0.5*90^2/30^2))./(2.*rsp + (rp + rn).*(1 + exp(-0.5*180^2/30^2)));

load my_fig3_trial5_fineRp.mat
%%
% centralPercentile2 prompt
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
fraction_rp1 = count_rp_total / experiment_num;
fraction_alpha1 = count_alpha_total / experiment_num;
fraction_rsp1 = count_rsp_total / experiment_num;
fraction_angle1 = count_angle_total / experiment_num;
fraction_sigma1 = count_sigma_total / experiment_num;
%%
writematrix(fraction_rp1,'Fig3_trials_5_fineRp.xlsx','Sheet','rp')
writematrix(fraction_alpha1,'Fig3_trials_5_fineRp.xlsx','Sheet','alpha')
writematrix(fraction_rsp1,'Fig3_trials_5_fineRp.xlsx','Sheet','rsp')
writematrix(fraction_angle1,'Fig3_trials_5_fineRp.xlsx','Sheet','angle')
writematrix(fraction_sigma1,'Fig3_trials_5_fineRp.xlsx','Sheet','sigma')
%%
% mode and central of 25%, 50% and 75%
offset = [0.25 0.5 0.75];
count_rp_total = zeros(numel(offset),cell_type);
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
        pdf_orientation_angle = pdf_angle(1:numel(pdf_angle)/2) + pdf_angle((numel(pdf_angle)/2 + 1):end); % combine directional angles to orientational angles
        pdf_orientation_angle(end+1) = pdf_orientation_angle(1);
        pdf_sigma = output(idx).marginal_likelihood.sigma.likelihoods;

        grid_rp = output(idx).marginal_likelihood.Rp.values;
        grid_alpha = output(idx).marginal_likelihood.Alpha.values;
        grid_rsp = output(idx).marginal_likelihood.Rsp.values;
        grid_angle = output(idx).marginal_likelihood.theta_pref.values;
        grid_orientation_angle = output(idx).marginal_likelihood.theta_pref.values(1:numel(pdf_orientation_angle));
        grid_sigma = output(idx).marginal_likelihood.sigma.values;

        [rp_bound,count_rp] = analysis.in_boundaries(grid_rp,pdf_rp,rp(i),offset); % Define offsets = 0.25, 0.5 and 0.75 of the central intervals
        [alpha_bound,count_alpha] = analysis.in_boundaries(grid_alpha,pdf_alpha,alpha(i),offset);
        [rsp_bound,count_rsp] = analysis.in_boundaries(grid_rsp,pdf_rsp,rsp(i),offset);
        [angle_bound,count_angle] = analysis.in_boundaries(grid_orientation_angle,pdf_orientation_angle,angle,offset);
        [sigma_bound,count_sigma] = analysis.in_boundaries(grid_sigma,pdf_sigma,sigma,offset);

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

    % plot cdf example curve
    % figure(11),
    % plot(unique_rp,unique_cdf_rp,'ro')
    % hold on,
    % plot(lower_rp,lower_bound_rp,'bo',upper_rp,upper_bound_rp,'bo'),
    % plot(cdf_rp_total,0:0.01:1,'k.')
    % hold off;
    % 
    % figure(12),
    % plot(unique_alpha,unique_cdf_alpha,'ro')
    % hold on,
    % plot(lower_alpha,lower_bound_alpha,'bo',upper_alpha,upper_bound_alpha,'bo'),
    % plot(cdf_alpha_total,0:0.01:1,'k.')
    % hold off;
    % 
    % figure(13),
    % plot(unique_rsp,unique_cdf_rsp,'ro')
    % hold on,
    % plot(lower_rsp,lower_bound_rsp,'bo',upper_rsp,upper_bound_rsp,'bo'),
    % plot(cdf_rsp_total,0:0.01:1,'k.')
    % hold off;
    % 
    % figure(14),
    % plot(unique_angle,unique_cdf_angle,'ro')
    % hold on,
    % plot(lower_angle,lower_bound_angle,'bo',upper_angle,upper_bound_angle,'bo'),
    % plot(cdf_angle_total,0:0.01:1,'k.')
    % hold off;
    % 
    % figure(15),
    % plot(unique_sigma,unique_cdf_sigma,'ro')
    % hold on,
    % plot(lower_sigma,lower_bound_sigma,'bo',upper_sigma,upper_bound_sigma,'bo'),
    % plot(cdf_sigma_total,0:0.01:1,'k.')
    % hold off;
%%
curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
plot(data_ideal{i}.angles,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',1)
errorbar(data_noisy{i}(58).angles,data_noisy{i}(58).mean_responses,data_noisy{i}(58).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*100-42).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',1)
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
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
figure(6),hold on;
for i = 1:5
plot(center,output(i*100-50).descriptors.oi.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
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
plot(output(i*100-50).marginal_likelihood.Rp.values,output(i*100-50).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
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
plot(output(i*100-50).marginal_likelihood.sigma.values,output(i*100-50).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
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
plot(output(i*100-50).marginal_likelihood.theta_pref.values,output(i*100-50).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
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
plot(output(i*100-50).marginal_likelihood.Rsp.values,output(i*100-50).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xticks(0:2:10),
yticks(0:.2:1),
set(gca,'FontSize',8)
xlabel('C','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('J','FontSize',12);

interval_tag = ['25%';'50%';'75%'];

figure(11),
hold on,
b = bar(offset,count_alpha_total,'grouped');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
xlabel('Credential Intervals of \alpha');
ylabel('Fraction of True Value in the Intervals (%)');
title('K','FontSize',12);

figure(12),
hold on,
b = bar(offset,count_rp_total,'grouped');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
xlabel('Credential Intervals of R_{p}');
ylabel('Fraction of True Value in the Intervals (%)');
title('L','FontSize',12);

figure(13),
hold on,
b = bar(offset,count_sigma_total,'grouped');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
xlabel('Credential Intervals of \sigma');
ylabel('Fraction of True Value in the Intervals (%)');
title('M','FontSize',12);

figure(14),
hold on,
b = bar(offset,count_angle_total,'grouped');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
xlabel('Credential Intervals of \theta_{pref}');
ylabel('Fraction of True Value in the Intervals (%)');
title('N','FontSize',12);

figure(15),
hold on,
b = bar(offset,count_rsp_total,'grouped');
for i = 1:5
 b(i).FaceColor = curvecolor{i};
end
xticks(offset);
xticklabels(interval_tag);
ylim([0,100]);
yticks(0:20:100);
xlabel('Credential Intervals of C');
ylabel('Fraction of True Value in the Intervals (%)');
title('O','FontSize',12);
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

path = 'D:\GitHub\vhlab-bayesoridir-matlab\Bayes_Estimation_paper\figures\noise_mdl\kcs\';

exportgraphics(F3a,[path 'Figure_3a_mode.pdf'],"ContentType","vector"),
exportgraphics(F3b,[path 'Figure_3b_mode.pdf'],"ContentType","vector"),
exportgraphics(F3c,[path 'Figure_3c_mode.pdf'],"ContentType","vector"),
exportgraphics(F3d,[path 'Figure_3d_mode.pdf'],"ContentType","vector"),
exportgraphics(F3e,[path 'Figure_3e_mode.pdf'],"ContentType","vector"),
exportgraphics(F3f,[path 'Figure_3f_mode.pdf'],"ContentType","vector"),
exportgraphics(F3g,[path 'Figure_3g_mode.pdf'],"ContentType","vector"),
exportgraphics(F3h,[path 'Figure_3h_mode.pdf'],"ContentType","vector"),
exportgraphics(F3i,[path 'Figure_3i_mode.pdf'],"ContentType","vector"),
exportgraphics(F3j,[path 'Figure_3j_mode.pdf'],"ContentType","vector"),
exportgraphics(F3k,[path 'Figure_3k_mode.pdf'],"ContentType","vector"),
exportgraphics(F3l,[path 'Figure_3l_mode.pdf'],"ContentType","vector"),
exportgraphics(F3m,[path 'Figure_3m_mode.pdf'],"ContentType","vector"),
exportgraphics(F3n,[path 'Figure_3n_mode.pdf'],"ContentType","vector"),
exportgraphics(F3o,[path 'Figure_3o_mode.pdf'],"ContentType","vector");