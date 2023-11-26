clear all,clc,close all;
%test5 : heatmap plotting for marginal likelihood of parameters
% extract file name
file = dir('*my_output_0*');
filename = {file.name};
%simulated 'true' curves
cell_num = 5;
rp = linspace(0,10,cell_num);
rn = linspace(0,5,cell_num);
rsp = linspace(10,0,cell_num);
% create plotting edges(n+1) and center(n) N = curves number
edges = 0:0.05:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
% calculate the 'true OI value' of all 50 simulate curve.
true_oi = 1 - (2.*rsp + 2.*(rp + rn).*exp(-0.5*90^2/30^2))./(2.*rsp + (rp + rn).*(1 + exp(-0.5*90^2/30^2)));

plot.get_rpref(true_oi,filename,cell_num,1)
plot.get_alpha(true_oi,filename,cell_num,1)
plot.get_rsp(true_oi,filename,cell_num,1)
plot.get_sigma(true_oi,filename,cell_num,1)
plot.get_thetapref(true_oi,filename,cell_num,1)