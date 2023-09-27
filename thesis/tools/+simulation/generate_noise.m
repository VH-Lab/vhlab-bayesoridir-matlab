function [noisy_curves,noise] = generate_noise(raw_curve,curve_num)
% generate <curve_num> tunning curves with random gaussian noise

%   Inputs: RAW_CURVE - 
%           CURVE_NUM - 
%
%   Output: NOISY_CURVE -
%           NOISE - gaussian distribution N(raw_datapoint,1)

if size(raw_curve,1) < size(raw_curve,2)
    raw_curve = raw_curve';
end
noisy_curves = zeros(length(raw_curve),curve_num);

noise = 0.5 .* raw_curve .* randn(size(noisy_curves));
noisy_curves = mean(raw_curve + noise);