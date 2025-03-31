clear,clc,close;
% load my_fig4_new_kcs.mat
% test6 -  prediction of different DI value
%simulated 'true' curves
% simulate sampling
cell_num = 5;
measurement_num = 5;
experiment_num = 100;
ang = 0:45:359;
%curves parameters
rp = 10;
rn = linspace(10,0,cell_num);
rsp = 0;
angle = 45;
sigma = 30;
%%
% bayes grid input
I = struct('Rp',linspace(0.1,20,41), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,31), ...
    'Sig',linspace(1,60,41), ...
    'Rsp',linspace(-10,10,41));
%%
%   generate simulate data
%   store in cell matrix. one cell one curve.
for i = 1:cell_num
    parameters1 = data.generate_fixed_parameters(rp,rn(i),angle,sigma,rsp);
    data_ideal{i} = data.generate_simulate_data(ang,parameters1);
    [data_mean{i}] = data.generate_noise(data_ideal{i},measurement_num,experiment_num);
    %plotting
    figure(),hold on,
    for j = 1:experiment_num
        plot(data_ideal{i}.angles,data_ideal{i}.responses,'k','LineWidth',2)
        plot(data_mean{i}(j).angles,data_mean{i}(j).mean_responses,'b*')
        ylim([0,15])
        xlabel('theta')
        ylabel('response magnitude')
        title('simulate tuning curve (+50% noise)')
    end
end
%%
%noise fitting model
m = [];
v = [];
for i = 1:cell_num
    for j = 1:experiment_num
        v = [v;data_mean{i}(j).responses_stddev];
        m = [m,data_mean{i}(j).mean_responses];
    end
end
figure(),
noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);

%%
% bayes fitting
tic
for i = 1:cell_num
    for j = 1:experiment_num
        fprintf('the fitting is at %d loop.\n',(i-1)*experiment_num + j)
        data_mean{i}(j).num_trials = measurement_num; % number of simulate trials
        [output((i-1)*experiment_num + j),~] = bayes_grid_function_proportional_noise_gpu(I,data_mean{i}(j),noise_mdl);
        toc;
    end
end
