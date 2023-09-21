clc; clear all;
load direction_sf_tf_data.mat
clear orientation_direction_bayes

cell_num = length(direction_data_s);
curve_num = 5;
% input structure of parameters
I = struct('Rp',logspace(log10(.5),log10(150),100), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,50), ...
    'Sig',logspace(log10(5),log10(90),20));
% fitting noise model
m = [];
v = [];
for ind = 1:cell_num
    v = [v,direction_data_s{ind}.document_properties.stimulus_tuningcurve.response_stddev];
    m = [m,direction_data_s{ind}.document_properties.stimulus_tuningcurve.response_mean];
end
index = (find(v>0));
v2 = v(index);
index = (find(m>0));
m2 = m(index);
figure(1),plot(log10(m2),log10(v2),'k.'),
xlabel('log10(mean)');
ylabel('log10(stddev)');
mdl = fitlm(log10(m2),log10(v2)),
figure(2),plot(mdl);
noise_coefficients = mdl.Coefficients{:,1};