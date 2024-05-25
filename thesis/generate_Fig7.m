clear all;close;clc;
%To access the records from Steve's old data, follow these examples:

%To check which cells exhibit strong responses before motion exposure, use:

load stevesolddata.mat

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


% To extract direction tuning data for each cell before motion exposure, use:

data_before = [];
figure(1),hold on;
for i=1:numel(good_indexes_before),
    for j=1:numel(cell_structures{good_indexes_before(i)}),
        if strcmp(cell_structures{good_indexes_before(i)}(j).type,'TP Ach OT Response curve'),
            data_before{end+1} = cell_structures{good_indexes_before(i)}(j).data;
            noisy_before(i).angle = data_before{end}(1,:);
            noisy_before(i).mean_responses = data_before{end}(2,:)';
            noisy_before(i).responses_stddev = data_before{end}(3,:);
            noisy_before(i).responses_stderr = data_before{end}(4,:);
            plot(noisy_before(i).angle,noisy_before(i).mean_responses)
        end;
    end;
end;
% To extract direction tuning data for each cell after motion exposure, use:

data_motionexposure = [];
figure(2),hold on;
for i=1:numel(good_indexes_motionexposure),
    for j=1:numel(cell_structures{good_indexes_motionexposure(i)}),
        if strcmp(cell_structures{good_indexes_motionexposure(i)}(j).type,'TP ME Ach OT Response curve'),
            data_motionexposure{end+1} = cell_structures{good_indexes_motionexposure(i)}(j).data;
            noisy_motionexposure(i).angle = data_motionexposure{end}(1,:);
            noisy_motionexposure(i).mean_responses = data_motionexposure{end}(2,:)';
            noisy_motionexposure(i).responses_stddev = data_motionexposure{end}(3,:);
            noisy_motionexposure(i).responses_stderr = data_motionexposure{end}(4,:);
            plot(noisy_motionexposure(i).angle,noisy_motionexposure(i).mean_responses,'.'),
        end;
    end;
end;

% To access the bootstrap data before motion exposure,

data_bootsbefore = [];
for i=1:numel(good_indexes_before),
    for j=1:numel(cell_structures{good_indexes_before(i)}),
        if strcmp(cell_structures{good_indexes_before(i)}(j).type,'TP Ach OT Bootstrap Carandini Fit Params'),
            data_bootsbefore{end+1} = cell_structures{good_indexes_before(i)}(j).data;
        end;
    end;
end;

for i = 1:numel(data_bootsbefore),
    rsp_min(i) = min(data_bootsbefore{i}(:,1));
    rsp_max(i) = max(data_bootsbefore{i}(:,1));
    rp_min(i) = min(data_bootsbefore{i}(:,2));
    rp_max(i) = max(data_bootsbefore{i}(:,2));
    op_min(i) = min(data_bootsbefore{i}(:,3));
    op_max(i) = max(data_bootsbefore{i}(:,3));
    sig_min(i) = min(data_bootsbefore{i}(:,4));
    sig_max(i) = max(data_bootsbefore{i}(:,4));
    rn_min(i) = min(data_bootsbefore{i}(:,5));
    rn_max(i) = max(data_bootsbefore{i}(:,5));
end;
rsp_range = [min(rsp_min),max(rsp_max)]
rp_range = [min(rp_min),max(rp_max)]
op_range = [min(op_min),max(op_max)]
sig_range = [min(sig_min),max(sig_max)]
rn_range = [min(rn_min),max(rn_max)]
% To access the bootstrap data after motion exposure,

data_bootsmotionexposure = [];
for i=1:numel(good_indexes_motionexposure),
    for j=1:numel(cell_structures{good_indexes_motionexposure(i)}),
        if strcmp(cell_structures{good_indexes_motionexposure(i)}(j).type,'TP ME Ach OT Bootstrap Carandini Fit Params'),
            data_bootsmotionexposure{end+1} = cell_structures{good_indexes_motionexposure(i)}(j).data;
        end;
    end;
end;

%%
% Bayesian Estimation
% bayes grid input
I = struct('Rp',linspace(0.001,0.3,30), ...
    'Op',0:30:359, ...
    'Alpha',linspace(0,1,11), ...
    'Sig',linspace(1,60,30), ...
    'Rsp',linspace(-0.1,0.2,20));
%%
%noise fitting model
m = [];
v = [];
for i = 1:numel(noisy_before)
    v = [v,noisy_before(i).responses_stderr];
    m = [m,noisy_before(i).mean_responses];
end
v = v(m>0);
m = m(m>0);
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
    [output(i),Lik] = bayes_grid_function_proportional_noise(I,noisy_before(i),noise_coefficients);
end
toc;


