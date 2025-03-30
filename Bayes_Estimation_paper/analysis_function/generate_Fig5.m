clear;clc;clf;
% generate figure 5
% load my_fig5_1_kcs_numtrials.mat

% simulate sampling
cell_num = 5;
measurement_num = 5;
experiment_num = 100;
ang{1} = 0:90:359;
ang{2} = 0:45:359;
ang{3} = 0:30:359;
ang{4} = 0:22.5:359;
ang{5} = 0:18:359;
%curves parameters
rp = 7.5;
rn = 3.75;
rsp = 2.5;
angle = 45;
sigma = 30;
%%
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%%
%   generate simulate data
%   store in cell matrix. one cell one curve.
for i = 1:cell_num
    parameters1 = data.generate_fixed_parameters(rp,rn,angle,sigma,rsp);
    data_ideal{i} = data.generate_simulate_data(ang{i},parameters1);
    [data_mean{i}] = data.generate_noise(data_ideal{i},measurement_num,experiment_num);
    %plotting
    figure(),hold on,
    for j = 1:experiment_num
        plot(data_ideal{i}.angle,data_ideal{i}.responses,'k','LineWidth',2)
        plot(data_mean{i}(j).angle,data_mean{i}(j).mean_responses,'b*')
        ylim([0,15])
        xlabel('theta')
        ylabel('response magnitude')
        title('simulate tuning curve (+50% noise)')
    end
end
%%
% Noise least-square fitting
m = [];
v = [];
for i = 1:numel(data_mean)
    for j = 1:numel(data_mean{i})
        m = [m;data_mean{i}(j).mean_responses(:)];
        v = [v;data_mean{i}(j).responses_stddev(:)];
    end
end

noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);

%%
% Bayes estimation
output = [];
for i = 1:numel(data_mean)
    for j = 1:numel(data_mean{i})
        data_mean{i}(j).num_trials = measurement_num; % number of simulate trials
        [output{end+1},~] = bayes_grid_function_proportional_noise_gpu(I,data_mean{i}(j),noise_mdl);
    end
end
clear i j;