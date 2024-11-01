clear;close all;clc;
% DATE（MMDDYYYY） - 06092024, Figure 8
% compare the result of Bayes estimation from neuron after motion
% exposure(AME) and before motion exposure(BME) and pick 3 typical cells to
% evaluate both Bayes Estimation and Bootstrap fitting results of
% theta_pref.

% Bayes fitting - 06082024
load my_fig8_kcs_BME_numtrials.mat

output_BME = output;

clear output I noisy_before

% Bayes fitting - 06082024
load my_fig8_kcs_AME_moreOp.mat

output_AME = output;

clear output I noisy_ME

load("stevesolddata.mat")

%%
% extract bootstrap result

good_indexes_before = [];
good_indexes_motionexposure = [];

for i=1:numel(cell_structures),
    for j=1:numel(cell_structures{i})
        if strcmp(cell_structures{i}(j).type,'TP Ach OT vec varies p')
            if cell_structures{i}(j).data<0.05
                good_indexes_before(end+1) = i;
            end
        end
        if strcmp(cell_structures{i}(j).type,'TP ME Ach OT vec varies p')
            if cell_structures{i}(j).data<0.05
                good_indexes_motionexposure(end+1) = i;
            end
        end
    end
end

good_both = intersect(good_indexes_before,good_indexes_motionexposure);
% data column 1: Response Offset
% column 2 : Rpref
% column 3 : theta_pref
% column 4 : sigma (tuning width)
% column 5 : Rnull
data_bootsBME = [];
data_bootsAME = [];
for i=1:numel(good_both)
    for j=1:numel(cell_structures{good_both(i)}),
        if strcmp(cell_structures{good_both(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params')
            data_bootsBME{i} = cell_structures{good_both(i)}(j).data;
        end
        if strcmp(cell_structures{good_both(i)}(j).type,'TP ME Ach OT Bootstrap Carandini Fit Params')
            data_bootsAME{i} = cell_structures{good_both(i)}(j).data;
        end
    end
end

%%
tic;
rsp_bootsBME = zeros(numel(data_bootsBME),1);
rp_bootsBME = zeros(numel(data_bootsBME),1);
op_bootsBME = zeros(numel(data_bootsBME),1);
sig_bootsBME = zeros(numel(data_bootsBME),1);
rn_bootsBME = zeros(numel(data_bootsBME),1);
rsp_BayesBME = zeros(numel(data_bootsBME),1);
rp_BayesBME = zeros(numel(data_bootsBME),1);
op_BayesBME = zeros(numel(data_bootsBME),1);
sig_BayesBME = zeros(numel(data_bootsBME),1);
rn_BayesBME = zeros(numel(data_bootsBME),1);

rsp_bootsAME = zeros(numel(data_bootsAME),1);
rp_bootsAME = zeros(numel(data_bootsAME),1);
op_bootsAME = zeros(numel(data_bootsAME),1);
sig_bootsAME = zeros(numel(data_bootsAME),1);
rn_bootsAME = zeros(numel(data_bootsAME),1);
rsp_BayesAME = zeros(numel(data_bootsAME),1);
rp_BayesAME = zeros(numel(data_bootsAME),1);
op_BayesAME = zeros(numel(data_bootsAME),1);
sig_BayesAME = zeros(numel(data_bootsAME),1);
rn_BayesAME = zeros(numel(data_bootsAME),1);

better_pick_cells = [6 9 10 27 30 32 34 35 36 37 39 40 43 49 52 54 59 61 66 67 70 77 78 85 91 96 97 101 105 107 120 124 130 143 151 155 156 166 167 171 182 186 188 189 195 197 207 209 211 215 217 220 226 229 248 262];
best_pick_cells = [120 40 32];
% Comparison of BME
for i = 1:numel(data_bootsBME),
    rsp_min(i) = min(data_bootsBME{i}(:,1));
    rsp_max(i) = max(data_bootsBME{i}(:,1));
    rp_min(i) = min(data_bootsBME{i}(:,2));
    rp_max(i) = max(data_bootsBME{i}(:,2));
    op_min(i) = min(data_bootsBME{i}(:,3));
    op_max(i) = max(data_bootsBME{i}(:,3));
    sig_min(i) = min(data_bootsBME{i}(:,4));
    sig_max(i) = max(data_bootsBME{i}(:,4));
    rn_min(i) = min(data_bootsBME{i}(:,5));
    rn_max(i) = max(data_bootsBME{i}(:,5));


    % bayes grid input
    var = max(noisy_BME(i).mean_responses);

    I_BME = struct('Rp',linspace(0.001,3*var,60), ...
        'Op',0:5:359, ...
        'Alpha',linspace(0,1,21), ...
        'Sig',linspace(1,60,60), ...
        'Rsp',linspace(-var,var,60), ...
        'Rn',linspace(0,3*var,60));

    rsp_N = histc(data_bootsBME{i}(:,1),linspace(-var,var,40));
    rsp_N = rsp_N./sum(rsp_N,'all');
    rp_N = histc(data_bootsBME{i}(:,2),linspace(0.001,3*var,60));
    rp_N = rp_N./sum(rp_N,'all');
    op_N = histc(data_bootsBME{i}(:,3),0:5:359);
    op_N = op_N./sum(op_N,'all');
    sig_N = histc(data_bootsBME{i}(:,4),linspace(1,60,60));
    sig_N = sig_N./sum(sig_N,'all');
    rn_N = histc(data_bootsBME{i}(:,5),linspace(0,3*var,60));
    rn_N = rn_N./sum(rn_N,'all');

    if ismember(i,best_pick_cells)
        figure(i),
        tiledlayout("vertical"),
        nexttile,hold on,
        plot(0:5:359,output_BME(i).marginal_likelihood.theta_pref.likelihoods,'b','LineWidth',1)
        plot(0:5:359,op_N,'b--','LineWidth',1)
        xlim([-5 365])
        ylim([0,1])
        set(gca,'FontSize',8)
        legend('Bayes BME','Bootstrap BME','FontSize',7)
        xlabel('\theta_{pref}','FontSize',10)
        ylabel('Likelihood','FontSize',10)
    end

    [Rp_value,Alpha_value] = meshgrid(I_BME.Rp,I_BME.Alpha);
    Rn_value = Alpha_value .* Rp_value;
    [Rn_value,ind_rn] = sort(Rn_value(:));

    [Rp_lik,Alpha_lik] = meshgrid(output_BME(i).marginal_likelihood.Rp.likelihoods,output_BME(i).marginal_likelihood.Alpha.likelihoods);
    Rn_lik = Rp_lik .* Alpha_lik;
    Rn_lik = Rn_lik(ind_rn);
    % Extract the indices of the Rn vector within specified interval ranges
    [~,ind_v] = histc(Rn_value,linspace(0,3*var,60));
    % Use the indices extracted from the numerical vector to accumulate the
    % probabilities within the same intervals in the probability vector
    summed_Rn_lik = accumarray(ind_v,Rn_lik);

    % Examine the likelihood distribution of Bayesian estimation and bootstrap with different parameters.
    % % bar(output_BME(i).marginal_likelihood.theta_pref.values,op_N);
    % % hold on,
    % % plot(output_BME(i).marginal_likelihood.theta_pref.values,output_BME(i).marginal_likelihood.theta_pref.likelihoods,'LineWidth',1),
    % % hold off;
    % % pause(1),

    % Obtain the peak values in the probability distribution of each parameter

    rsp_BayesBME(i) = I_BME.Rsp(output_BME(i).marginal_likelihood.Rsp.likelihoods == max(output_BME(i).marginal_likelihood.Rsp.likelihoods));
    rp_BayesBME(i) = I_BME.Rp(output_BME(i).marginal_likelihood.Rp.likelihoods == max(output_BME(i).marginal_likelihood.Rp.likelihoods));
    op_BayesBME(i) = I_BME.Op(output_BME(i).marginal_likelihood.theta_pref.likelihoods == max(output_BME(i).marginal_likelihood.theta_pref.likelihoods));
    sig_BayesBME(i) = I_BME.Sig(output_BME(i).marginal_likelihood.sigma.likelihoods == max(output_BME(i).marginal_likelihood.sigma.likelihoods));
    rn_BayesBME(i) = I_BME.Rn(summed_Rn_lik == max(summed_Rn_lik));

if numel(I_BME.Rsp(rsp_N == max(rsp_N))) == 1
        rsp_bootsBME(i) = I_BME.Rsp(rsp_N == max(rsp_N));
    else
        rsp_bootsBME(i) = median(data_bootsBME{i}(:,1)); % current solution: taking median value
    end

    if numel(I_BME.Rp(rp_N == max(rp_N))) == 1
        rp_bootsBME(i) = I_BME.Rp(rp_N == max(rp_N));
    else
        rp_bootsBME(i) = median(data_bootsBME{i}(:,2));
    end

    if numel(I_BME.Op(op_N == max(op_N))) == 1
        op_bootsBME(i) = I_BME.Op(op_N == max(op_N));
    else
        op_bootsBME(i) = median(data_bootsBME{i}(:,3)); 
    end

    if numel(I_BME.Sig(sig_N == max(sig_N))) == 1
        sig_bootsBME(i) = I_BME.Sig(sig_N == max(sig_N));
    else
        sig_bootsBME(i) = median(data_bootsBME{i}(:,4));
    end

    if numel(I_BME.Rn(rn_N == max(rn_N))) == 1
        rn_bootsBME(i) = I_BME.Rn(rn_N == max(rn_N));
    else
        rn_bootsBME(i) = median(data_bootsBME{i}(:,5));
    end
end


% comparison of AME
for i = 1:numel(data_bootsAME),
    rsp_min(i) = min(data_bootsAME{i}(:,1));
    rsp_max(i) = max(data_bootsAME{i}(:,1));
    rp_min(i) = min(data_bootsAME{i}(:,2));
    rp_max(i) = max(data_bootsAME{i}(:,2));
    op_min(i) = min(data_bootsAME{i}(:,3));
    op_max(i) = max(data_bootsAME{i}(:,3));
    sig_min(i) = min(data_bootsAME{i}(:,4));
    sig_max(i) = max(data_bootsAME{i}(:,4));
    rn_min(i) = min(data_bootsAME{i}(:,5));
    rn_max(i) = max(data_bootsAME{i}(:,5));


    % bayes grid input
    var = max(noisy_AME(i).mean_responses);

    I_AME = struct('Rp',linspace(0.001,3*var,60), ...
        'Op',0:1:359, ...
        'Alpha',linspace(0,1,21), ...
        'Sig',linspace(1,60,60), ...
        'Rsp',linspace(-var,var,60), ...
        'Rn',linspace(0,3*var,60));

    rsp_N = histc(data_bootsAME{i}(:,1),linspace(-var,var,40));
    rsp_N = rsp_N./sum(rsp_N,'all');
    rp_N = histc(data_bootsAME{i}(:,2),linspace(0.001,3*var,60));
    rp_N = rp_N./sum(rp_N,'all');
    op_N = histc(data_bootsAME{i}(:,3),0:1:359);
    op_N = op_N./sum(op_N,'all');
    sig_N = histc(data_bootsAME{i}(:,4),linspace(1,60,60));
    sig_N = sig_N./sum(sig_N,'all');
    rn_N = histc(data_bootsAME{i}(:,5),linspace(0,3*var,60));
    rn_N = rn_N./sum(rn_N,'all');

    if ismember(i,best_pick_cells)
        figure(i),
        nexttile,hold on,
        plot(0:1:359,output_AME(i).marginal_likelihood.theta_pref.likelihoods,'r','LineWidth',1)
        plot(0:1:359,op_N,'r--','LineWidth',1)
        xlim([-5 365])
        ylim([0,1])
        set(gca,'FontSize',8)
        legend('Bayes AME','Bootstrap AME','FontSize',7)
        xlabel('\theta_{pref}','FontSize',10)
        ylabel('Likelihood','FontSize',10)
    end

    [Rp_value,Alpha_value] = meshgrid(I_AME.Rp,I_AME.Alpha);
    Rn_value = Alpha_value .* Rp_value;
    [Rn_value,ind_rn] = sort(Rn_value(:));

    [Rp_lik,Alpha_lik] = meshgrid(output_AME(i).marginal_likelihood.Rp.likelihoods,output_AME(i).marginal_likelihood.Alpha.likelihoods);
    Rn_lik = Rp_lik .* Alpha_lik;
    Rn_lik = Rn_lik(ind_rn);
    % Extract the indices of the Rn vector within specified interval ranges
    [~,ind_v] = histc(Rn_value,linspace(0,3*var,60));
    % Use the indices extracted from the numerical vector to accumulate the
    % probabilities within the same intervals in the probability vector
    summed_Rn_lik = accumarray(ind_v,Rn_lik);

    % Examine the likelihood distribution of Bayesian estimation and bootstrap with different parameters.
    % % bar(output_BME(i).marginal_likelihood.theta_pref.values,op_N);
    % % hold on,
    % % plot(output_BME(i).marginal_likelihood.theta_pref.values,output_BME(i).marginal_likelihood.theta_pref.likelihoods,'LineWidth',1),
    % % hold off;
    % % pause(1),

    % Obtain the peak values in the probability distribution of each parameter

    rsp_BayesAME(i) = I_AME.Rsp(output_AME(i).marginal_likelihood.Rsp.likelihoods == max(output_AME(i).marginal_likelihood.Rsp.likelihoods));
    rp_BayesAME(i) = I_AME.Rp(output_AME(i).marginal_likelihood.Rp.likelihoods == max(output_AME(i).marginal_likelihood.Rp.likelihoods));
    op_BayesAME(i) = I_AME.Op(output_AME(i).marginal_likelihood.theta_pref.likelihoods == max(output_AME(i).marginal_likelihood.theta_pref.likelihoods));
    sig_BayesAME(i) = I_AME.Sig(output_AME(i).marginal_likelihood.sigma.likelihoods == max(output_AME(i).marginal_likelihood.sigma.likelihoods));
    rn_BayesAME(i) = I_AME.Rn(summed_Rn_lik == max(summed_Rn_lik));

    if numel(I_AME.Rsp(rsp_N == max(rsp_N))) == 1
        rsp_bootsAME(i) = I_AME.Rsp(rsp_N == max(rsp_N));
    else
        rsp_bootsAME(i) = median(data_bootsAME{i}(:,1)); % current solution: taking mean value
    end

    if numel(I_AME.Rp(rp_N == max(rp_N))) == 1
        rp_bootsAME(i) = I_AME.Rp(rp_N == max(rp_N));
    else
        rp_bootsAME(i) = median(data_bootsAME{i}(:,2));
    end

    if numel(I_AME.Op(op_N == max(op_N))) == 1
        op_bootsAME(i) = I_AME.Op(op_N == max(op_N));
    else
        op_bootsAME(i) = median(data_bootsAME{i}(:,3));
    end

    if numel(I_AME.Sig(sig_N == max(sig_N))) == 1
        sig_bootsAME(i) = I_AME.Sig(sig_N == max(sig_N));
    else
        sig_bootsAME(i) = median(data_bootsAME{i}(:,4));
    end

    if numel(I_AME.Rn(rn_N == max(rn_N))) == 1
        rn_bootsAME(i) = I_AME.Rn(rn_N == max(rn_N));
    else
        rn_bootsAME(i) = median(data_bootsAME{i}(:,5));
    end
end
toc;
%%
figure(1),hold on,
rsp_min_axis = -.24;
rsp_max_axis = .21;
plot([rsp_min_axis,rsp_max_axis],[rsp_min_axis,rsp_max_axis],'k--')
plot(rsp_bootsBME,rsp_BayesBME,'sb','LineWidth',1)
plot(rsp_bootsAME,rsp_BayesAME,'xr','LineWidth',1)
xlim([rsp_min_axis,rsp_max_axis])
ylim([rsp_min_axis,rsp_max_axis])
xticks(-.2:.1:.21)
yticks(-.2:.1:.21)
set(gca,'FontSize',8)
legend('Before Motion Exposure','After Motion Exposure','Location','northwest','FontSize',7)
xlabel('Mean of C from Bootstrap','FontSize',10)
ylabel('ML of C from Bayes','FontSize',10)
title('G','FontSize',12);

figure(2),hold on,
rp_min_axis = -0.01;
rp_max_axis = 0.6;
plot([rp_min_axis,rp_max_axis],[rp_min_axis,rp_max_axis],'k--')
plot(rp_bootsBME,rp_BayesBME,'sb','LineWidth',1)
plot(rp_bootsAME,rp_BayesAME,'xr','LineWidth',1)
xlim([rp_min_axis,rp_max_axis])
ylim([rp_min_axis,rp_max_axis])
xticks(0:.1:.6)
yticks(0:.1:.6)
set(gca,'FontSize',8)
xlabel('Mean of R_{p} from Bootstrap','FontSize',10)
ylabel('ML of R_{p} from Bayes','FontSize',10)
title('H','FontSize',12);

figure(3),hold on,
op_min_axis = -5;
op_max_axis = 365;
plot([op_min_axis,op_max_axis],[op_min_axis,op_max_axis],'k--')
plot(op_bootsBME,op_BayesBME,'sb','LineWidth',1)
plot(op_bootsAME,op_BayesAME,'xr','LineWidth',1)
xlim([op_min_axis,op_max_axis])
ylim([op_min_axis,op_max_axis])
xticks(0:100:400)
yticks(0:100:400)
set(gca,'FontSize',8)
xlabel('Mean of \theta_{pref} from Bootstrap','FontSize',10)
ylabel('ML of \theta_{pref} from Bayes','FontSize',10)
title('I','FontSize',12);

figure(4),hold on,
sig_min_axis = -1;
sig_max_axis = 62;
plot([sig_min_axis,sig_max_axis],[sig_min_axis,sig_max_axis],'k--')
plot(sig_bootsBME,sig_BayesBME,'sb','LineWidth',1)
plot(sig_bootsAME,sig_BayesAME,'xr','LineWidth',1)
xlim([sig_min_axis,sig_max_axis])
ylim([sig_min_axis,sig_max_axis])
xticks(0:10:60)
yticks(0:10:60)
set(gca,'FontSize',8)
xlabel('Mean of \sigma from Bootstrap','FontSize',10)
ylabel('ML of \sigma from Bayes','FontSize',10)
title('J','FontSize',12);

figure(5),hold on,
rn_min_axis = -.01;
rn_max_axis = .5;
plot([rn_min_axis,rn_max_axis],[rn_min_axis,rn_max_axis],'k--')
plot(rn_bootsBME,rn_BayesBME,'sb','LineWidth',1)
plot(rn_bootsAME,rn_BayesAME,'xr','LineWidth',1)
xlim([rn_min_axis,rn_max_axis])
ylim([rn_min_axis,rn_max_axis])
xticks(0:.1:.5)
yticks(0:.1:.5)
set(gca,'FontSize',8)
legend('Before Motion Exposure','After Motion Exposure','Location','southeast','FontSize',7)
xlabel('Mean of R_{n} from Bootstrap','FontSize',10)
ylabel('ML of R_{n} from Bayes','FontSize',10)
title('K','FontSize',12);
%%
% export figures
% F8a = figure(32);
% F8b = figure(40);
% F8c = figure(120);
% F8d = figure(1);
% F8e = figure(2);
% F8f = figure(3);
% F8g = figure(4);
% F8h = figure(5);
% 
% path = 'D:\GitHub\vhlab-bayesoridir-matlab\thesis\figures\noise_mdl\kcs\';
% 
% exportgraphics(F8a,[path 'Figure_8ab_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8b,[path 'Figure_8cd_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8c,[path 'Figure_8ef_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8d,[path 'Figure_8g_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8e,[path 'Figure_8h_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8f,[path 'Figure_8i_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8g,[path 'Figure_8j_kcs_numtrials.pdf'],"ContentType","vector"),
% exportgraphics(F8h,[path 'Figure_8k_kcs_numtrials.pdf'],"ContentType","vector"),
