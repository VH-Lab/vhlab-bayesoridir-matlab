clear all,clc,clf;

%% DATE（MMDDYYYY） - 06092024
% compare the result of Bayes estimation from neuron before motion exposure(BME)

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
% data column 1: Response Offset
% column 2 : Rpref
% column 3 : theta_pref
% column 4 : sigma (tuning width)
% column 5 : Rnull
data_bootsBME = [];
data_bootsAME = [];
for i=1:numel(good_both),
    for j=1:numel(cell_structures{good_both(i)}),
        if strcmp(cell_structures{good_both(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params'),
            data_bootsBME{i} = cell_structures{good_both(i)}(j).data;
        end;
        if strcmp(cell_structures{good_both(i)}(j).type,'TP ME Ach OT Bootstrap Carandini Fit Params'),
            data_bootsAME{i} = cell_structures{good_both(i)}(j).data;
        end;
    end;
end;

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
% comparison of BME
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
    elseif ismember(rsp_BayesBME(i),I_BME.Rsp(rsp_N == max(rsp_N)))
        rsp_bootsBME(i) = rsp_BayesBME(i);
    else
        rsp_bootsBME(i) = mean(I_BME.Rsp(rsp_N == max(rsp_N))); % current solution: taking mean value
    end

    if numel(I_BME.Rp(rp_N == max(rp_N))) == 1
        rp_bootsBME(i) = I_BME.Rp(rp_N == max(rp_N));
    elseif ismember(rp_BayesBME(i),I_BME.Rp(rp_N == max(rp_N)))
        rp_bootsBME(i) = rp_BayesBME(i);
    else
        rp_bootsBME(i) = mean(I_BME.Rp(rp_N == max(rp_N)));
    end

    if numel(I_BME.Op(op_N == max(op_N))) == 1
        op_bootsBME(i) = I_BME.Op(op_N == max(op_N));
    elseif ismember(op_BayesBME(i),I_BME.Op(op_N == max(op_N)))
        op_bootsBME(i) = op_BayesBME(i);
    else
        op_bootsBME(i) = mean(I_BME.Op(op_N == max(op_N)));
        I_BME.Op(op_N == max(op_N))
    end

    if numel(I_BME.Sig(sig_N == max(sig_N))) == 1
        sig_bootsBME(i) = I_BME.Sig(sig_N == max(sig_N));
    elseif ismember(sig_BayesBME(i),I_BME.Sig(sig_N == max(sig_N)))
        sig_bootsBME(i) = sig_BayesBME(i);
    else
        sig_bootsBME(i) = mean(I_BME.Sig(sig_N == max(sig_N)));
    end

    if numel(I_BME.Rn(rn_N == max(rn_N))) == 1
        rn_bootsBME(i) = I_BME.Rn(rn_N == max(rn_N));
    elseif ismember(rn_BayesBME(i),I_BME.Rn(rn_N == max(rn_N)))
        rn_bootsBME(i) = rn_BayesBME(i);
    else
        rn_bootsBME(i) = mean(I_BME.Rn(rn_N == max(rn_N)));
    end
end
toc;
%%
figure(1),
plot(rsp_bootsBME,rsp_BayesBME,'*')
rsp_min_axis = 1.05 * min([rsp_bootsBME(:);rsp_BayesBME(:)]);
rsp_max_axis = 1.05 * max([rsp_bootsBME(:);rsp_BayesBME(:)]);
xlim([rsp_min_axis,rsp_max_axis])
ylim([rsp_min_axis,rsp_max_axis])
xlabel('Rsp ETA from Bootstrap')
ylabel('Rsp ETA from Bayes Estimation')

figure(2),
plot(rp_bootsBME,rp_BayesBME,'*')
rp_min_axis = 1.05 * min([rp_bootsBME(:);rp_BayesBME(:)]);
rp_max_axis = 1.05 * max([rp_bootsBME(:);rp_BayesBME(:)]);
xlim([rp_min_axis,rp_max_axis])
ylim([rp_min_axis,rp_max_axis])
xlabel('Rp ETA from Bootstrap')
ylabel('Rp ETA from Bayes Estimation')

figure(3),
plot(op_bootsBME,op_BayesBME,'*')
op_min_axis = min([op_bootsBME(:);op_BayesBME(:)]) - 0.05 * max([op_bootsBME(:);op_BayesBME(:)]);
op_max_axis = 1.05 * max([op_bootsBME(:);op_BayesBME(:)]);
xlim([op_min_axis,op_max_axis])
ylim([op_min_axis,op_max_axis])
xlabel('Op ETA from Bootstrap')
ylabel('Op ETA from Bayes Estimation')

figure(4),
plot(sig_bootsBME,sig_BayesBME,'*')
sig_min_axis = 1.05 * min([sig_bootsBME(:);sig_BayesBME(:)]);
sig_max_axis = 1.05 * max([sig_bootsBME(:);sig_BayesBME(:)]);
xlim([sig_min_axis,sig_max_axis])
ylim([sig_min_axis,sig_max_axis])
xlabel('Sigma ETA from Bootstrap')
ylabel('Sigma ETA from Bayes Estimation')

figure(5),
plot(rn_bootsBME,rn_BayesBME,'*')
rn_min_axis = min([rn_bootsBME(:);rn_BayesBME(:)]) - 0.05 * max([rn_bootsBME(:);rn_BayesBME(:)]);
rn_max_axis = 1.05 * max([rn_bootsBME(:);rn_BayesBME(:)]);
xlim([rn_min_axis,rn_max_axis])
ylim([rn_min_axis,rn_max_axis])
xlabel('Rn ETA from Bootstrap')
ylabel('Rn ETA from Bayes Estimation')
