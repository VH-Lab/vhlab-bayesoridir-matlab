clear;close;clc;
% DATE（MMDDYYYY） - 06082024
% Making Bayes estimation for neurons before and after motion exposure 

%extract experiment data
load stevesolddata.mat,
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
% To extract direction tuning data and number of trials for each cell before motion exposure, use:
data_BME = [];
for i=1:numel(good_both),
    for j=1:numel(cell_structures{good_both(i)}),
        if strcmp(cell_structures{good_both(i)}(j).type,'TP Ach OT Response struct'),
           data_trials_num = numel(cell_structures{good_both(i)}(j).data.ind{1});% gain number of trials at each direction
           noisy_BME(i).num_trials = data_trials_num;
        end
        if strcmp(cell_structures{good_both(i)}(j).type,'TP Ach OT Response curve'),
            data_BME{end+1} = cell_structures{good_both(i)}(j).data;
            noisy_BME(i).angle = data_BME{end}(1,:);
            noisy_BME(i).mean_responses = data_BME{end}(2,:)';
            noisy_BME(i).responses_stddev = data_BME{end}(3,:);
            noisy_BME(i).responses_stderr = data_BME{end}(4,:);
        end;
    end;
end;

%noise fitting model
m = [];
v = [];
for i = 1:numel(noisy_BME),
    v = [v;noisy_BME(i).responses_stderr(:)];
    m = [m;noisy_BME(i).mean_responses(:)];
end

noise_mdl = vis.bayes.noise.fit_proportional_noise_plus_c(m,v,1);

v1 = v(m>0);
m1 = m(m>0);
mdl = fitlm(log10(m1),log10(v1)),
figure(),plot(mdl);
xlabel('log10(response mean)'),
ylabel('log10(response stddev)'),
noise_coefficients = mdl.Coefficients{:,1};



%%
% bayes fitting
tic;
for i = 1:length(noisy_BME),
    % bayes grid input
var = max(noisy_BME(i).mean_responses);

I = struct('Rp',linspace(0.001,3*var,60), ...
    'Op',0:1:359, ...
    'Alpha',linspace(0,1,21), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(-var,var,60));

    fprintf('the fitting is at %d loop of %d.\n',i,numel(noisy_BME)),
    [output(i),~] = bayes_grid_function_proportional_noise_gpu(I,noisy_BME(i),noise_mdl);
end
toc;