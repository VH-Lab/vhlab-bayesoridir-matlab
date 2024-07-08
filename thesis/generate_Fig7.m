clear all;close;clc;
%To access the records from Steve's old data, follow these examples:

%To check which cells exhibit strong responses before motion exposure, use:

load stevesolddata.mat
good_indexes_BME = [];
good_indexes_AME = [];

for i=1:numel(cell_structures),
    for j=1:numel(cell_structures{i}),
        if strcmp(cell_structures{i}(j).type,'TP Ach OT vec varies p'),
            if cell_structures{i}(j).data<0.05,
                good_indexes_BME(end+1) = i;
            end;
        end;
        if strcmp(cell_structures{i}(j).type,'TP ME Ach OT vec varies p'),
            if cell_structures{i}(j).data<0.05,
                good_indexes_AME(end+1) = i;
            end;
        end;
    end;
end;


% To extract direction tuning data for each cell before motion exposure, use:

data_BME = [];
noisy_BME = [];
for i=1:numel(good_indexes_BME),
    for j=1:numel(cell_structures{good_indexes_BME(i)}),
        if strcmp(cell_structures{good_indexes_BME(i)}(j).type,'TP Ach OT Response curve'),
            data_BME{end+1} = cell_structures{good_indexes_BME(i)}(j).data;
            noisy_BME(i).angle = data_BME{end}(1,:);
            noisy_BME(i).mean_responses = data_BME{end}(2,:)';
            noisy_BME(i).responses_stddev = data_BME{end}(3,:);
            noisy_BME(i).responses_stderr = data_BME{end}(4,:);
        end;
    end;
end;
% To extract direction tuning data for each cell after motion exposure, use:

data_AME = [];
noisy_AME = [];
for i=1:numel(good_indexes_AME),
    for j=1:numel(cell_structures{good_indexes_AME(i)}),
        if strcmp(cell_structures{good_indexes_AME(i)}(j).type,'TP ME Ach OT Response curve'),
            data_AME{end+1} = cell_structures{good_indexes_AME(i)}(j).data;
            noisy_AME(i).angle = data_AME{end}(1,:);
            noisy_AME(i).mean_responses = data_AME{end}(2,:)';
            noisy_AME(i).responses_stddev = data_AME{end}(3,:);
            noisy_AME(i).responses_stderr = data_AME{end}(4,:);
        end;
    end;
end;

% To access the bootstrap data before motion exposure,

data_bootsBME = [];
for i=1:numel(good_indexes_BME),
    for j=1:numel(cell_structures{good_indexes_BME(i)}),
        if strcmp(cell_structures{good_indexes_BME(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params'),
            data_bootsBME{end+1} = cell_structures{good_indexes_BME(i)}(j).data;
        end;
    end;
end;

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
end;
rsp_range = [min(rsp_min),max(rsp_max)]
rp_range = [min(rp_min),max(rp_max)]
op_range = [min(op_min),max(op_max)]
sig_range = [min(sig_min),max(sig_max)]
rn_range = [min(rn_min),max(rn_max)]
% To access the bootstrap data after motion exposure,

data_bootsAME = [];
for i=1:numel(good_indexes_AME),
    for j=1:numel(cell_structures{good_indexes_AME(i)}),
        if strcmp(cell_structures{good_indexes_AME(i)}(j).type,'TP ME Ach OT Bootstrap Carandini Fit Params'),
            data_bootsAME{end+1} = cell_structures{good_indexes_AME(i)}(j).data;
        end;
    end;
end;

% Bayesian Estimation
% bayes grid input
I = struct('Rp',linspace(0.001,0.3,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,21), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(-0.1,0.2,60));

%noise fitting model
m = [];
v = [];
for i = 1:numel(noisy_BME)
    v = [v,noisy_BME(i).responses_stderr];
    m = [m,noisy_BME(i).mean_responses];
end

v = v(m > 0.01);
m = m(m > 0.01);
mdl = fitlm(log10(m),log10(v)),
figure(),plot(mdl);
xlabel('log10(response mean)'),
ylabel('log10(response stddev)')
noise_coefficients = mdl.Coefficients{:,1};
%%
% bayes fitting
tic;
for i = 1:5,
    fprintf('the fitting is at %d loop.\n',i),
    [output(i),~] = bayes_grid_function_proportional_noise_gpu(I,noisy_BME(i),noise_coefficients);
end
toc;


