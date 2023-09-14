function test()
%TEST tests the acurracy and performance of vis.bayes.double_gaussian.grid_proportional_noise
%   Detailed explanation goes here


c.Rsp = [ 1 ];
c.Rp = 10;
c.Rn = 3;
c.theta_pref = 45;
c.sigma = 30;
c.noise_offset = 0.4;
c.noise_slope = 0.4;

c.data.angles = [ 0:30:360-30]';
c.data.num_trials = 5*ones(size(c.data.angles));
c.data.mean_responses_raw = vlt.math.double_gauss_180(c.data.angles,[c.Rsp c.Rp c.theta_pref c.sigma c.Rn]);
c.data.noise_sigma = vis.bayes.noise.proportional([c.noise_offset c.noise_slope],c.data.mean_responses_raw,c.data.num_trials);
c.data.mean_responses = c.data.mean_responses_raw + randn(size(c.data.mean_responses_raw)) .* c.data.noise_sigma;

grid = struct('Rp',logspace(log10(.5),log10(150),100), ...
    'Op',0:5:359, ...
    'Alpha',linspace(0,1,50), ...
    'Sig',logspace(log10(5),log10(90),10), ...
    'Rsp',logspace(log10(0.01),log10(100),5));


[output_struct,lik] = vis.bayes.double_gaussian.grid_proportional_noise(grid,c.data,[c.noise_offset c.noise_slope]);

