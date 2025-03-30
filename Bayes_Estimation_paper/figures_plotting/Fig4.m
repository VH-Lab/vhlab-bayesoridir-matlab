clear;close;clc;
load my_fig4_mode.mat
% Fig 4 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning(DI).

cell_num = 5;
experiment_num = 100;
ang = 0:45:359;
%curves parameters
rp = 10;
rn = linspace(10,0,cell_num);
alpha = rn./rp;
rsp = 0;
angle = 45;
sigma = 30;
%%
% mode and central of 25%, 50% and 75%
offset = [0.25 0.5 0.75];
count_rp_total = zeros(numel(offset),cell_num);
count_alpha_total = count_rp_total;
count_rsp_total = count_rp_total;
count_angle_total = count_rp_total;
count_sigma_total = count_rp_total;

% Get all pdf in output data.
for i = 1:cell_num
    for j = 1:experiment_num
        idx = (i-1)*experiment_num + j;
        pdf_rp = output(idx).marginal_likelihood.Rp.likelihoods;
        pdf_alpha = output(idx).marginal_likelihood.Alpha.likelihoods;
        pdf_rsp = output(idx).marginal_likelihood.Rsp.likelihoods;
        pdf_angle = output(idx).marginal_likelihood.theta_pref.likelihoods;
        pdf_sigma = output(idx).marginal_likelihood.sigma.likelihoods;

        grid_rp = output(j).marginal_likelihood.Rp.values;
        grid_alpha = output(j).marginal_likelihood.Alpha.values;
        grid_rsp = output(j).marginal_likelihood.Rsp.values;
        grid_angle = output(j).marginal_likelihood.theta_pref.values;
        grid_sigma = output(j).marginal_likelihood.sigma.values;

        [rp_bound,count_rp] = analysis.in_boundaries(grid_rp,pdf_rp,rp,offset); % Define offsets = 0.25, 0.5 and 0.75 of the central intervals
        [alpha_bound,count_alpha] = analysis.in_boundaries(grid_alpha,pdf_alpha,alpha(i),offset);
        [rsp_bound,count_rsp] = analysis.in_boundaries(grid_rsp,pdf_rsp,rsp,offset);
        [angle_bound,count_angle] = analysis.in_boundaries(grid_angle,pdf_angle,angle,offset);
        [sigma_bound,count_sigma] = analysis.in_boundaries(grid_sigma,pdf_sigma,sigma,offset);

        count_rp_total(:,i) = count_rp_total(:,i) + count_rp(:);
        count_alpha_total(:,i) = count_alpha_total(:,i) + count_alpha(:);
        count_rsp_total(:,i) = count_rsp_total(:,i) + count_rsp(:);
        count_angle_total(:,i) = count_angle_total(:,i) + count_angle(:);
        count_sigma_total(:,i) = count_sigma_total(:,i) + count_sigma(:);
        % Calculate alpha interval boundaries
        % Include the boundary effect when the lower or upper interval boundary
        % hit the edge of cdf.

        % Calculate rsp interval boundaries
        % Include the boundary effect when the lower or upper interval boundary
        % hit the edge of cdf.

        % Calculate angle interval boundaries
        % Include the boundary effect when the lower or upper interval boundary
        % hit the edge of cdf.

        % Calculate sigma interval boundaries
        % Include the boundary effect when the lower or upper interval boundary
        % hit the edge of cdf.
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
% calculate the 'true DI value' from all 50 simulated curves.
true_di = (rp - rn).*(1 - exp(-0.5*180^2./sigma^2))./(rsp + rp + rn.*exp(-0.5.*180^2./sigma^2));
% extract di likelihood data in output structure. Transform the data form into
% plotting
for i = 1:5
parameters1 = data.generate_fixed_parameters(rp,rn(i),angle,sigma,rsp);
end

curvecolor = {'#0072BD','#77AC30','#EDB120','#7E2F8E','#D95319'};
titlename = {'A','B','C','D','E'};
for i = 1:5
figure(i),hold on;
plot(data_ideal{i}.angle,data_ideal{i}.responses,'Color',curvecolor{i},'LineWidth',1)
errorbar(data_mean{i}(10).angle,data_mean{i}(10).mean_responses,data_mean{i}(10).responses_stderr,'*','Color','#A2142F','MarkerSize',7)
plot(0:359,output(i*100-90).maximum_likelihood.parameters.tunning_curve,'k','LineWidth',1)
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
plot(center,output(i*10-3).descriptors.di.histogram_likelihoods,'Color',curvecolor{i},'LineWidth',1)
end
ylim([-0.02,1]),
xticks(0:.2:1),
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('DI','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('F','FontSize',12);

figure(7),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.Rp.values,output(i*10-3).marginal_likelihood.Rp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
end
ylim([-0.02,1]),
xticks(0:5:20), 
yticks(0:.2:1),
set(gca,'FontSize',8),
xlabel('R_{p}','FontSize',10),
ylabel('Likelihood','FontSize',10),
title('G','FontSize',12);

figure(8),
hold on,
for i = 1:5
plot(output(i*10-3).marginal_likelihood.sigma.values,output(i*10-3).marginal_likelihood.sigma.likelihoods,'Color',curvecolor{i},'LineWidth',1),
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
plot(output(i*10-3).marginal_likelihood.theta_pref.values,output(i*10-3).marginal_likelihood.theta_pref.likelihoods,'Color',curvecolor{i},'LineWidth',1),
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
plot(output(i*10-3).marginal_likelihood.Rsp.values,output(i*10-3).marginal_likelihood.Rsp.likelihoods,'Color',curvecolor{i},'LineWidth',1),
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
% F4a = figure(1);
% F4b = figure(2);
% F4c = figure(3);
% F4d = figure(4);
% F4e = figure(5);
% F4f = figure(6);
% F4g = figure(7);
% F4h = figure(8);
% F4i = figure(9);
% F4j = figure(10);
% 
% path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs\';
% 
% exportgraphics(F4a,[path 'Figure_4a_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4b,[path 'Figure_4b_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4c,[path 'Figure_4c_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4d,[path 'Figure_4d_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4e,[path 'Figure_4e_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4f,[path 'Figure_4f_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4g,[path 'Figure_4g_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4h,[path 'Figure_4h_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4i,[path 'Figure_4i_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F4j,[path 'Figure_4j_kcs_numtrials.pdf'],"ContentType","vector"),