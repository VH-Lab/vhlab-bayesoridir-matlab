clear all,clc,close all;
% test1 - prediction of different OI value
%simulated 'true' various curves / single curve
cell_type = 5;
measurement_num = 5;
experiment_num = 10;
ang = 0:45:359;
rp = linspace(0,10,5);
rn = linspace(0,5,5);
rsp = linspace(10,0,5);
%%
% bayes grid input
I = struct('Rp',linspace(0.1,20,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,15), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(0.1,10,60));
%%
%   generate simulate data
%   store in cell matrix. one cell one curve.
for i = 1:cell_type
    parameters1 = data.generate_fixed_parameters(rp(i),rn(i),45,30,rsp(i));
    data_ideal{i} = data.generate_simulate_data(ang,parameters1);
    [data_noisy{i}] = data.generate_noise(data_ideal{i},measurement_num,experiment_num);
    %plotting
    figure(),hold on,
    for j = 1:experiment_num
        plot(data_ideal{i}.angle,data_ideal{i}.responses,'k','LineWidth',2)
        errorbar(data_noisy{i}(j).angle,data_noisy{i}(j).mean_responses,data_noisy{i}(j).responses_stderr,'b*')
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
for i = 1:cell_type
    for j = 1:experiment_num
        v = [v;data_noisy{i}(j).responses_stddev];
        m = [m,data_noisy{i}(j).mean_responses];
    end
end
mdl = fitlm(log10(m),log10(v)),
figure(),plot(mdl);
xlabel('log10(response mean)'),
ylabel('log10(response stddev)')  
noise_coefficients = mdl.Coefficients{:,1};
%%
% bayes fitting
tic
for i = 1:cell_type
    for j = 1:experiment_num
        fprintf('the fitting is at %d loop.\n',(i-1)*experiment_num + j)
        [output((i-1)*experiment_num + j),~] = bayes_grid_function_proportional_noise_gpu(I,data_noisy{i}(j),noise_coefficients);
        toc;
    end
end
 
