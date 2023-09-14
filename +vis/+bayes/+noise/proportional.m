function [sigma] = proportional(noise_model, responses, num_trials)
%PROPORTIONAL computes proportional noise model
% 
% SIGMA = vis.bayes.noise.proportional(NOISE_MODEL, RESPONSES, NUM_TRIALS)
%
% Computes the standard deviations of a noise process that depends on
% the response in a manner such that it is linear on a log-log plot.
% 
% For a given response r and number of trials t, the sigma s is
%    s = (10^offset*r^slope)/sqrt(t)
%
% That is, it is assumed that the response is the average of the individual
% response to a number of presentations (num_trials) of a stimulus. The
% log-log fit is assumed to be generated from the standard deviation of
% individual trials (on the Y axis). Therefore, by the central limit
% theorum, the standard deviation of the average is the classic standard
% error of the mean (standard deviation / sqrt(t) ).
%
% Inputs:
%    NOISE_MODEL = [ offset slope]
%        where offset and slope are the parameters of a line in a LOG/LOG
%        plot of the standard deviation of the response (Y axis) against
%        the response itself (X axis).
%    RESPONSES - an N-dimensional vector of averaged response values to use for each
%        value of sigma to be created
%    NUM_TRIALS - an N-dimensional vector the number of independent trials
%        that were averaged to produce RESPONSES
% Output:
%    SIGMA - an N-dimensional vector of sigma values calculated according
%            to the formula above

if ~vlt.data.eqlen(size(responses),size(num_trials)),
    error(['responses and num_trials must have same dimensions.']);
end;


offset = noise_model(1);
slope = noise_model(2);

sigma = (10.^offset + responses.^slope)./sqrt(num_trials);

