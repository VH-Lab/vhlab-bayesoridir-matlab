clear all,clc,close all;
%test5 : heatmap plotting for marginal likelihood of parameters
% extract file name
file = dir('*my_output_0*');
filename = {file.name};
%simulated 'true' curves
cell_type = 5;
rp = 10;
rn = linspace(10,0,cell_type);
rsp = 0;
sigma = 30;
% create plotting edges(n+1) and center(n) N = curves number
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
% calculate the 'true OI value' of all 50 simulate curve.
true_di = (rp-rn) * (1-exp(-0.5*180^2./sigma^2))./(rsp + rp + rn*exp(-0.5*180^2./sigma^2));


plot.get_rpref(true_di,filename,cell_type,0)
plot.get_alpha(true_di,filename,cell_type,0)
plot.get_rsp(true_di,filename,cell_type,0)
plot.get_sigma(true_di,filename,cell_type,0)
plot.get_thetapref(true_di,filename,cell_type,0)