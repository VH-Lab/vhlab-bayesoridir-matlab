function [noisy_data] = generate_von_mises_noise(raw_curve,measurement_num,experiment_num)
% generate <curve_num> tunning curves with random gaussian noise
%
%   Inputs: RAW_CURVE - 'true' simulated tuning curve
%           MEASUREMENT_NUM - Repeat <measurement_num> simulated measurements in one experiment
%           EXPERIMENT_NUM - Repeat the same experiment <experiment_num> times
%   
%   Output: NOISY_DATA - N(0,1)

if size(raw_curve.responses,1) < size(raw_curve.responses,2)
    raw_curve.responses = raw_curve.responses(:);
end

noisy_responses = zeros(length(raw_curve.responses),experiment_num,measurement_num);
noise = normrnd(2,1,size(noisy_responses));
noisy_responses = raw_curve.responses + noise;
noisy_responses(noisy_responses<0) = 0; % Limitation: Response frequency must larger or equal to 0
noisy_curve.mean_responses = mean(noisy_responses,3);
noisy_curve.responses_stddev = std(noisy_responses,[],3);
noisy_curve.responses_stderr = std(noisy_responses,[],3)./repmat(sqrt(measurement_num),length(raw_curve.responses),experiment_num);

for i = 1:experiment_num
    noisy_data(i).angle = raw_curve.angle;
    noisy_data(i).mean_responses = noisy_curve.mean_responses(:,i)';
    noisy_data(i).responses_stddev = noisy_curve.responses_stddev(:,i);
    noisy_data(i).responses_stderr = noisy_curve.responses_stderr(:,i);
end

end