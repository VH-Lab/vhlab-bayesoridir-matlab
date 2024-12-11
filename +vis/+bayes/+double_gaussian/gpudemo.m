function [output_struct,lik,debug] = gpudemo(mode)

rsp_values = sort([-logspace(log10(0.1),log10(40),20) logspace(log10(0.1),log10(40),20)]);
rp_values = logspace(log10(0.1),log10(150),100);
alpha_values = 0:0.05:1;
thetap_values = 0:2:359;
sig_values = 10:5:90;

noise_model = [1.3549 1.9182 0.5461];

angles = [ 0:30:360-30 ];
P = [ 0.5 100 0 45 30 ];
resps = vis.oridir.doublegaussianfunc(angles,P);

param_grid = struct('Rsp',rsp_values,...
    'Rp',rp_values,...
    'Alpha',alpha_values,...
    'Op',thetap_values,...
    'Sig',sig_values);

resp_struct = struct('angles',angles(:),...
    'mean_responses',resps(:),...
    'num_trials',5*ones(size(resps(:))));

switch mode,
	case 'gpu-full',
		[output_struct,lik,debug] = vis.bayes.double_gaussian.grid_proportional_noise_metal_full(param_grid, resp_struct, noise_model);
	case 'gpu-fast',
		[output_struct,lik,debug] = vis.bayes.double_gaussian.grid_proportional_noise_metal(param_grid, resp_struct, noise_model);
	case 'gpu-pc',
		debug = [];
		[output_struct,lik] = vis.bayes.double_gaussian.grid_proportional_noise_gpu(param_grid, resp_struct, noise_model);
	case 'nogpu',
		debug = [];
		[output_struct,lik] = vis.bayes.double_gaussian.grid_proportional_noise(param_grid, resp_struct, noise_model,'useCPU',1,'verbose',1);
	otherwise,
		error(['Unknown method ' mode '.']);
end;

