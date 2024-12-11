# Bayesian parameter estimation of orientation and direction data

As described in Wu and Van Hooser (2024/2025), this library performs a full Bayesian parameter estimation of the parameters of a double gaussian fit. This is achieved by evaluating a range of possible parameter values. This operation is relatively slow on CPUs, but is fast on GPUs, and we offer calls to Apple Metal and NVIDIA GPUs on Linux/Windows systems.

## Model

The parameters that are estimated are $Rsp$, $Rp$, $Rn$, $\theta_{pref}$, and $\sigma$ in the following equation:

$$R(\theta) = Rsp + Rp * \text{exp}(-\text{angdiff}(\theta-\theta_{pref})^2/2\sigma^2) + Rn * \text{exp}(-\text{angdiff}(\theta-\theta_{pref}+180)^2/2\sigma^2)$$

where $\text{angdiff}(\theta)$ is the absolute angular difference, $Rsp$ is the non-specific response, $Rp$ is the peak response above the non-specific response to the preferred direction $\theta_{pref}$, $Rn$ is the response above the non-specific response to the opposite direction, and $\sigma$ is the width of the Gaussian.

## Requirements

- [vhlab-toolbox-matlab](https://github.com/VH-Lab/vhlab-toolbox-matlab)
- [MatlabMetal](https://github.com/stevevanhooser/MatlabMetal)
- Make sure that vhlab-bayesoridir-matlab is on the path

## Examples

1. Small example to explore the parameter estimation:

```[matlab]

% example equation
angles = [ 0:30:360-30 ];
P = [ 0.5 10 5 45 30 ]; % [Rsp Rp Rn theta_pref sigma] 
resps = vis.oridir.doublegaussianfunc(angles,P); % raw data

  %the parameter values to explore
rsp_values = sort([-logspace(log10(0.1),log10(40),20) logspace(log10(0.1),log10(40),20)]);
rp_values = logspace(log10(0.1),log10(150),20);
alpha_values = 0:0.05:1;
thetap_values = 0:5:359;
sig_values = 10:5:90;

% We use a noise model where the expected noise of the spike response
% depends on the mean response of the neuron
% noise model: N = C + K * M^S, where M is the mean response of neuron
noise_model = [1.3549 1.9182 0.5461]; % [C K S]

param_grid = struct('Rsp',rsp_values,...
    'Rp',rp_values,...
    'Alpha',alpha_values,...
    'Op',thetap_values,...
    'Sig',sig_values);

resp_struct = struct('angles',angles(:),...
    'mean_responses',resps(:),...
    'num_trials',5*ones(size(resps(:))));

[output_struct,lik] = vis.bayes.double_gaussian.grid_proportional_noise(param_grid, resp_struct, noise_model,'verbose',1);
  % or force running on CPU
%[output_struct,lik] = vis.bayes.double_gaussian.grid_proportional_noise(param_grid, resp_struct, noise_model,'useCPU',1,'verbose',1);

vis.bayes.double_guassian.plot_results(output_struct);
```

2. Research-grade parameter space exploration example (GPU strongly recommended):

```[matlab]

% example equation
angles = [ 0:30:360-30 ];
P = [ 0.5 10 5 45 30 ]; % [Rsp Rp Rn theta_pref sigma] 
resps = vis.oridir.doublegaussianfunc(angles,P); % raw data

  %the parameter values to explore
rsp_values = sort([-logspace(log10(0.1),log10(40),20) logspace(log10(0.1),log10(40),20)]);
rp_values = logspace(log10(0.1),log10(150),100);
alpha_values = 0:0.05:1;
thetap_values = 0:2:359;
sig_values = 10:5:90;

% We use a noise model where the expected noise of the spike response
% depends on the mean response of the neuron
% noise model: N = C + K * M^S, where M is the mean response of neuron
noise_model = [1.3549 1.9182 0.5461]; % [C K S]

param_grid = struct('Rsp',rsp_values,...
    'Rp',rp_values,...
    'Alpha',alpha_values,...
    'Op',thetap_values,...
    'Sig',sig_values);

resp_struct = struct('angles',angles(:),...
    'mean_responses',resps(:),...
    'num_trials',5*ones(size(resps(:))));

[output_struct,lik] = vis.bayes.double_gaussian.grid_proportional_noise(param_grid, resp_struct, noise_model,'verbose',1);

vis.bayes.double_guassian.plot_results(output_struct);
```

## Material related to the paper

The code to produce all figures of Wu and Van Hooser 2024/2025 is in the folder `Bayes_Estimation_paper`.


