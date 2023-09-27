function [noisy_curve,noise] = generate_noise(raw_curve,curve_num)
% generate <curve_num> tunning curves with random gaussian noise

%   Inputs: RAW_CURVE - 
%           CURVE_NUM - 
%
%   Output: NOISY_CURVE - 
%           NOISE - N(raw_datapoint,1)

if size(raw_curve,1) < size(raw_curve,2)
    raw_curve = raw_curve';
end
noisy_curve = zeros(length(raw_curve),curve_num);

noise = 0.5 .* raw_curve .* randn(size(noisy_curve));
noisy_curve = mean(raw_curve + noise,2);
end