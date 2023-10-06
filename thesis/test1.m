clear,clc;

%simulated 'true' curves
ang = 0:45:359;
rp = linspace(0,10,5);
rn = linspace(0,5,5);
rsp = linspace(10,0,5);
measurement_num = 5;
experiment_num = 10;

% bayes grid input
I = struct('Rp',logspace(log10(.5),log10(150),100), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,50), ...
    'Sig',logspace(log10(5),log10(90),10), ...
    'Rsp',logspace(log10(0.01),log10(100),5));
%%
%   generate simulate data
%   store in cell matrix. one cell one curve.
for i = 1:5
    parameters1 = data.generate_fixed_parameters(rp(i),rn(i),45,30,rsp(i));
    data_ideal{i} = data.generate_simulate_data(ang,parameters1);
    [data_mean{i}] = data.generate_noise(data_ideal{i},measurement_num,experiment_num);
    %plotting
    figure(),hold on,
    plot(data_ideal{i}.angle,data_ideal{i}.responses,'k','LineWidth',2)
    plot(data_mean{i}.angle,data_mean{i}.mean_responses,'b')
    ylim([0,15])
    xlabel('theta')
    ylabel('response magnitude')
    title('simulate tuning curve (+50% noise)')
end
%%
%noise fitting model
m = [];
v = [];
for i = 1:5
    v = [v;data_mean{i}.responses_stddev];
    m = [m,data_mean{i}.mean_responses];
end
mdl = fitlm(log10(m),log10(v)),
figure(),plot(mdl);
xlabel('log10(response mean)'),
ylabel('log10(response stddev)')  
noise_coefficients = mdl.Coefficients{:,1};
%%
% bayes fitting
tic;
for i = 1:5
[output(i),Lik] = bayes_grid_function_proportional_noise(I,data_mean{i},noise_coefficients);
end
toc;

