clear;clc,clf;
%% DATE（MMDDYYYY） - 06092024, Figure 8
% compare the result of Bayes estimation from neuron after motion
% exposure(AME) and before motion exposure(BME) and pick 3 typical cells to
% evaluate both Bayes Estimation and Bootstrap fitting results of
% theta_pref.

load("bayes_before_motionexposure.mat")

output_BME = output;
noisy_BME = noisy_before;
clear output I noisy_before

load("bayes_afterME.mat")

output_AME = output;
noisy_AME = noisy_ME;

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
        'Rsp',linspace(-var,var,40), ...
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
        plot(0:5:359,output_BME(i).marginal_likelihood.theta_pref.likelihoods,'b')
        plot(0:5:359,op_N,'b--')
        ylim([0,1])
        legend('Bayes BME','Bootstrap BME')
        xlabel('\theta_{pref}')
        ylabel('Likelihood')
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
        'Op',0:5:359, ...
        'Alpha',linspace(0,1,21), ...
        'Sig',linspace(1,60,60), ...
        'Rsp',linspace(-var,var,40), ...
        'Rn',linspace(0,3*var,60));

    rsp_N = histc(data_bootsAME{i}(:,1),linspace(-var,var,40));
    rsp_N = rsp_N./sum(rsp_N,'all');
    rp_N = histc(data_bootsAME{i}(:,2),linspace(0.001,3*var,60));
    rp_N = rp_N./sum(rp_N,'all');
    op_N = histc(data_bootsAME{i}(:,3),0:5:359);
    op_N = op_N./sum(op_N,'all');
    sig_N = histc(data_bootsAME{i}(:,4),linspace(1,60,60));
    sig_N = sig_N./sum(sig_N,'all');
    rn_N = histc(data_bootsAME{i}(:,5),linspace(0,3*var,60));
    rn_N = rn_N./sum(rn_N,'all');

    if ismember(i,best_pick_cells)
        figure(i),
        nexttile,hold on,
        plot(0:5:359,output_AME(i).marginal_likelihood.theta_pref.likelihoods,'r')
        plot(0:5:359,op_N,'r--')
        ylim([0,1])
        legend('Bayes AME','Bootstrap AME')
        xlabel('\theta_{pref}')
        ylabel('Likelihood')
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
plot(rsp_bootsBME,rsp_BayesBME,'sb')
plot(rsp_bootsAME,rsp_BayesAME,'xr')
rsp_min_axis = 1.05 * min([rsp_bootsAME(:);rsp_BayesAME(:);rsp_bootsBME(:);rsp_BayesBME(:)]);
rsp_max_axis = 1.05 * max([rsp_bootsAME(:);rsp_BayesAME(:);rsp_bootsBME(:);rsp_BayesBME(:)]);
plot([rsp_min_axis,rsp_max_axis],[rsp_min_axis,rsp_max_axis],'k--')
xlim([rsp_min_axis,rsp_max_axis])
ylim([rsp_min_axis,rsp_max_axis])
legend('Before Motion Exposure','After Motion Exposure','Location','northwest')
xlabel('Response Offset ETA From Bootstrap')
ylabel('Response Offset ETA From Bayes Estimation')

figure(2),hold on,
plot(rp_bootsBME,rp_BayesBME,'sb')
plot(rp_bootsAME,rp_BayesAME,'xr')
plot([0 360],[0 360],'k--')
rp_min_axis = 1.05 * min([rp_bootsAME(:);rp_BayesAME(:);rp_bootsBME(:);rp_BayesBME]);
rp_max_axis = 1.05 * max([rp_bootsAME(:);rp_BayesAME(:);rp_bootsBME(:);rp_BayesBME]);
plot([rp_min_axis,rp_max_axis],[rp_min_axis,rp_max_axis],'k--')
xlim([rp_min_axis,rp_max_axis])
ylim([rp_min_axis,rp_max_axis])
xlabel('R_{pref} ETA From Bootstrap')
ylabel('R_{pref} ETA From Bayes Estimation')

figure(3),hold on,
plot(op_bootsBME,op_BayesBME,'sb')
plot(op_bootsAME,op_BayesAME,'xr')
op_min_axis = min([op_bootsAME(:);op_BayesAME(:);op_bootsBME(:);op_BayesBME(:)]) - 0.05 * max([op_bootsAME(:);op_BayesAME(:);op_bootsBME(:);op_BayesBME(:)]);
op_max_axis = 1.05 * max([op_bootsAME(:);op_BayesAME(:);op_bootsBME(:);op_BayesBME(:)]);
plot([op_min_axis,op_max_axis],[op_min_axis,op_max_axis],'k--')
xlim([op_min_axis,op_max_axis])
ylim([op_min_axis,op_max_axis])
xlabel('\theta_{pref} ETA From Bootstrap')
ylabel('\theta_{pref} ETA from Bayes Estimation')

figure(4),hold on,
plot(sig_bootsBME,sig_BayesBME,'sb')
plot(sig_bootsAME,sig_BayesAME,'xr')
sig_min_axis = 1.05 * min([sig_bootsAME(:);sig_BayesAME(:);sig_bootsBME(:);sig_BayesBME(:)]);
sig_max_axis = 1.05 * max([sig_bootsAME(:);sig_BayesAME(:);sig_bootsBME(:);sig_BayesBME(:)]);
plot([sig_min_axis,sig_max_axis],[sig_min_axis,sig_max_axis],'k--')
xlim([sig_min_axis,sig_max_axis])
ylim([sig_min_axis,sig_max_axis])
xlabel('\sigma ETA From Bootstrap')
ylabel('\sigma ETA From Bayes Estimation')

figure(5),hold on,
plot(rn_bootsBME,rn_BayesBME,'sb')
plot(rn_bootsAME,rn_BayesAME,'xr')
rn_min_axis = min([rn_bootsAME(:);rn_BayesAME(:);rn_bootsBME(:);rn_BayesBME(:)]) - 0.05 * max([rn_bootsAME(:);rn_BayesAME(:);rn_bootsBME(:);rn_BayesBME(:)]);
rn_max_axis = 1.05 * max([rn_bootsAME(:);rn_BayesAME(:);rn_bootsBME(:);rn_BayesBME(:)]);
plot([rn_min_axis,rn_max_axis],[rn_min_axis,rn_max_axis],'k--')
xlim([rn_min_axis,rn_max_axis])
ylim([rn_min_axis,rn_max_axis])
xlabel('R_{null} ETA From Bootstrap')
ylabel('R_{null} ETA From Bayes Estimation')
