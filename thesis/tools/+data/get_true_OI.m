function [OI_underlying] = get_true_OI(parameters)

oi(rp,alpha,sig,rsp) = 1 - (2*parameters.response + 2*parameters.rpref *(1+parameters.alpha)*exp(-0.5*90^2/parameters.sigma^2))./(2*parameters.response + parameters.rpref *(1+parameters.alpha)*(1 + exp(-0.5*90^2/parameters.sigma ^2)));

OI_underlying = oi;

end