function [noisy_curve] = generate_noise(raw_curve,measurement_num,experiment_num)
% generate <curve_num> tunning curves with random gaussian noise
%
%   Inputs: RAW_CURVE - 'true' simulated tuning curve
%           MEASUREMENT_NUM - Repeat <measurement_num> simulated measurements in one experiment
%           EXPERIMENT_NUM - Repeat the same experiment <experiment_num> times
%   
%   Output: NOISY_CURVE - N(0,(0.5*rp)^2)

noisy_curve.angle = raw_curve.angle;

if size(raw_curve.responses,1) < size(raw_curve.responses,2)
    raw_curve.responses = raw_curve.responses';
end

noisy_responses = zeros(length(raw_curve.responses),experiment_num,measurement_num);
noise = 0.5 .* raw_curve.responses .* randn(size(noisy_responses));
noisy_responses = mean(raw_curve.responses + noise,3);

noisy_curve.mean_responses = mean(noisy_responses,2);
noisy_curve.responses_stddev = std(noisy_responses,[],2);
noisy_curve.responses_stderr = stderr(noisy_responses');

if size(noisy_curve.mean_responses,1) > size(noisy_curve.mean_responses,2)
    noisy_curve.mean_responses = noisy_curve.mean_responses';
end

end