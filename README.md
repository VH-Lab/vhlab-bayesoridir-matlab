# Bayesian parameter estimation of orientation and direction data

As described in Wu and Van Hooser (2024/2025), this library performs a full Bayesian parameter estimation of the parameters of a double gaussian fit. This is achieved by evaluating a range of possible parameter values. This operation is relatively slow on CPUs, but is fast on GPUs, and we offer calls to Apple Metal and NVIDIA GPUs on Linux/Windows systems.

## Model

The parameters that are estimated are $Rsp$, $Rp$, $Rn$, $\theta_{pref}$, and $\sigma$ in the following equation:

$$R(\theta) = Rsp + Rp * exp(-\text{angdiff}(\theta-\theta_{pref})^2/2\sigma^2) + Rn * exp(-\text{angdiff}(\theta-\theta_{pref}+180)^2/2\sigma^2)

## Requirements

- [vhlab-toolbox-matlab](https://github.com/VH-Lab/vhlab-toolbox-matlab)
- [MatlabMetal](https://github.com/stevevanhooser/MatlabMetal) 

## Examples

To run a small calculation, 
