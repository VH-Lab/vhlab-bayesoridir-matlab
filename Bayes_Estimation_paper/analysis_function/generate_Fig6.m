clear,clc,close all;

% generate figure 6
cell_num = 5;
measurement_num = 5;
experiment_num = 10;
ang = 0:45:359;
%curves parameters
rp = linspace(10,0,cell_num);
angle = 45;
sigma = 1;

%%
%   generate simulate data
%   store in cell matrix. one cell one curve.
for i = 1:cell_num
    parameters = data.vonmises.generate_von_mises_fixed_parameters(rp(i),angle,sigma);
    data_ideal{i} = data.vonmises.generate_von_mises_data(ang,parameters);
    [data_noisy{i}] = data.vonmises.generate_von_mises_noise(data_ideal{i},measurement_num,experiment_num);
    %plotting
    figure(),hold on,
    for j = 1:experiment_num
        plot(data_ideal{i}.angle,data_ideal{i}.responses,'k','LineWidth',2)
        plot(data_noisy{i}(j).angle,data_noisy{i}(j).mean_responses,'b*')
        ylim([0,15])
        xlabel('theta')
        ylabel('response magnitude')
        title('simulate tuning curve (+50% noise)')
    end
end
%%
m = [];
v = [];
for i = 1:cell_num
    for j = 1:experiment_num
        v = [v;data_noisy{i}(j).responses_stddev];
        m = [m,data_noisy{i}(j).mean_responses];
    end
end
noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);

% Bayes Estimate Grid Range and Size
I = struct('Rp',linspace(0.001,10,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.001,10,50));
%%
tic
for i = 1:cell_num
    for j = 1:experiment_num
        fprintf('the fitting is at %d loop.\n',(i-1)*experiment_num + j)
        data_noisy{i}(j).num_trials = measurement_num; % number of simulate trials
        [output((i-1)*experiment_num + j),~] = bayes_grid_function_proportional_noise_gpu(I,data_noisy{i}(j),noise_mdl);
        toc;
    end
end
