clear all;close;clc;
% DATE（MMDDYYYY） - 06082024
% Making Bayes estimation for neurons before and after motion exposure 

%extract experiment data
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

good_both = intersect(good_indexes_before,good_indexes_motionexposure);
% To extract direction tuning data for each cell before motion exposure, use:

data_ME = [];
for i=1:numel(good_both),
    for j=1:numel(cell_structures{good_both(i)}),
        if strcmp(cell_structures{good_both(i)}(j).type,'TP ME Ach OT Response curve'),
            data_ME{end+1} = cell_structures{good_both(i)}(j).data;
            noisy_ME(i).angle = data_ME{end}(1,:);
            noisy_ME(i).mean_responses = data_ME{end}(2,:)';
            noisy_ME(i).responses_stddev = data_ME{end}(3,:);
            noisy_ME(i).responses_stderr = data_ME{end}(4,:);
        end;
    end;
end;
% Bayesian Estimation

%noise fitting model
m = [];
v = [];
for i = 1:numel(noisy_ME)
    v = [v,noisy_ME(i).responses_stderr];
    m = [m,noisy_ME(i).mean_responses];
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
for i = 1:length(noisy_ME),
    % bayes grid input
var = max(noisy_ME(i).mean_responses);

I = struct('Rp',linspace(0.001,3*var,60), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,21), ...
    'Sig',linspace(1,60,60), ...
    'Rsp',linspace(-var,var,40));

    fprintf('the fitting is at %d loop.\n',i),
    [output(i),~] = bayes_grid_function_proportional_noise_gpu(I,noisy_ME(i),noise_coefficients);
end
toc;