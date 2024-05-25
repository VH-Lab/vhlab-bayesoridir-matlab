 clear all;close;clc;
load my_fig4.mat;
% Fig 4 - Performance of Bayesian parameter estimation for simulated data
% of varying orientation tuning.

cell_num = 5;
experiment_num = 10;
ang = 0:45:359;
%curves parameters
rp = 10;
rn = linspace(10,0,cell_num);
rsp = 0;
angle = 45;
sigma = 30;
% create plotting edges(n+1) and center(n) N = curves number
edges = 0:0.01:1;
center = edges(1:end-1) + (edges(1)+edges(2))/2;
% calculate the 'true DI value' from all 50 simulated curves.
true_di = (rp - rn).*(1 - exp(-0.5*180^2./sigma^2))./(rsp + rp + rn.*exp(-0.5.*180^2./sigma^2));             
% extract di likelihood data in output structure. Transform the data form into
% plotting suitable form
di_lik = zeros(length(output),length(output(1).descriptors.di.histogram_likelihoods));
for i = 1:length(output)
    di_lik(i,:) = output(i).descriptors.di.histogram_likelihoods;
end

di_lik_heatmap = zeros(length(center),length(true_di));
for i = 1:cell_num       
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_di_lik = mean(di_lik(ind_lower:ind_upper,:),1);
di_lik_heatmap(:,i) = avg_di_lik;
end
% extend y axie: A row of x
di_lik_pcolor = [di_lik_heatmap;zeros(1,length(true_di))];
%extend x axie: A column of y
di_lik_pcolor = [di_lik_pcolor,zeros(length(edges),1)];
%%
%plotting
% x-axie represent the value of 'true DI'
figure(1),
% tiledlayout flow
% nexttile
heatmap(true_di,flip(center),flip(di_lik_heatmap),'Colormap',parula)
xlabel('Value of ^\primetrue DI^\prime')
ylabel('Probability of DI Value')
title('Distribution of DI Probability')

% nexttile
% pcolor([true_di 1],edges,di_lik_pcolor)
% xlabel('Value of ^\primetrue DI^\prime')
% ylabel('Probability of DI Value')
% title('Distribution of DI Probability')
%%
% plotting by another way
% x_axie represent the bin of the value of 'true DI'
% put 'true DI' value into specific DI bin.
ind_di = discretize(true_di,edges);
% extract DI likelihood data in output structure. Each row represent one
% cell
di_lik = zeros(length(output),length(output(1).descriptors.di.histogram_likelihoods));
for i = 1:length(output)
    di_lik(i,:) = output(i).descriptors.di.histogram_likelihoods;
end
% creat a matrix to store the information of Heatmap grid
% predict DI value: 0:0.05:1
% likelihood of DI value: 0:0.05:1
di_lik_heatmap = zeros(length(output(1).descriptors.di.histogram_count),length(output(1).descriptors.di.histogram_likelihoods));
for i = 1:cell_num
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_di_lik = mean(di_lik(ind_lower:ind_upper,:),1);
% DI_lik_cumulative = discretize(avg_DI_lik,edges);
di_lik_heatmap(:,ind_di(i)) = avg_di_lik;
end
% extend y axie: should be the number of x (simulate 'true DI' bin)
di_lik_pcolor = [di_lik_heatmap;zeros(1,length(output(1).descriptors.di.histogram_count))];
%extend x axie: should be the number of y (likelihood bin)
di_lik_pcolor = [di_lik_pcolor,zeros(length(edges),1)];

figure(2),
% tiledlayout flow
% nexttile
heatmap(center,flip(center),flip(di_lik_heatmap),'Colormap',parula)
xlabel('Value of DI')
ylabel('Predict DI Value')
title('Distribution of DI Probability')

% nexttile
% pcolor(edges,edges,DI_lik_pcolor)
% xlabel('Value of DI')
% ylabel('Predict of DI Value')
% title('Distribution of DI Probability')
%%
% extract rpref likelihood data in output structure. Transform the data form into
% plotting suitable form
rp_lik = zeros(length(output),length(output(1).marginal_likelihood.Rp.likelihoods));
rp_value = output(1).marginal_likelihood.Rp.values;
for i = 1:length(output)
    rp_lik(i,:) = output(i).marginal_likelihood.Rp.likelihoods;
end

rp_lik_heatmap = zeros(length(rp_value),length(true_di));
for i = 1:cell_num
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_rp_lik = mean(rp_lik(ind_lower:ind_upper,:),1);
rp_lik_heatmap(:,i) = avg_rp_lik;
end
% extend y axie: should be the number of x (simulate 'true DI' bin)
rp_lik_pcolor = [rp_lik_heatmap;zeros(1,length(true_di))];
%extend x axie: should be the number of y (likelihood bin)
rp_lik_pcolor = [rp_lik_pcolor,zeros(1+length(rp_value),1)];
figure(3),
% tiledlayout flow
% nexttile
heatmap(true_di,flip(output(1).marginal_likelihood.Rp.values),flip(rp_lik_heatmap),'Colormap',parula)
xlabel('Value of ^\primetrue DI^\prime')
ylabel('Predict of Rpref Value')
title('Distribution of Rpref Probability')

% nexttile
% pcolor([true_di 1],[0 rp_value],rp_lik_pcolor)
% xlabel('Value of ^\primetrue DI^\prime')
% ylabel('Predict of Rpref Value')
% title('Distribution of Rpref Probability')

%%
% extract alpha likelihood data in output structure, making a plotting suitable form
alpha_lik = zeros(length(output),length(output(1).marginal_likelihood.Alpha.likelihoods));
alpha_value = output(1).marginal_likelihood.Alpha.values;
for i = 1:length(output)
    alpha_lik(i,:) = output(i).marginal_likelihood.Alpha.likelihoods;
end

alpha_lik_heatmap = zeros(length(alpha_value),length(true_di));
for i = 1:cell_num
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_alpha_lik = mean(alpha_lik(ind_lower:ind_upper,:),1);
alpha_lik_heatmap(:,i) = avg_alpha_lik;
end
% extend y axie: should be the number of x (simulate 'true DI' bin)
alpha_lik_pcolor = [alpha_lik_heatmap;zeros(1,length(true_di))];
%extend x axie: should be the number of y (likelihood bin)
alpha_lik_pcolor = [alpha_lik_pcolor,zeros(1+length(alpha_value),1)];
figure(4),
% tiledlayout flow
% nexttile
heatmap(true_di,flip(output(1).marginal_likelihood.Alpha.values),flip(alpha_lik_heatmap),'Colormap',parula)
xlabel('Value of ^\primetrue DI^\prime')
ylabel('Predict of Alpha Value')
title('Distribution of Alpha Probability')

% nexttile
% pcolor([true_di 1],[0 alpha_value],alpha_lik_pcolor)
% xlabel('Value of ^\primetrue DI^\prime')
% ylabel('Predict of Alpha Value')
% title('Distribution of Alpha Probability')
%%
% extract rsp likelihood data in output structure. Transform the data form into
% plotting suitable form
rsp_lik = zeros(length(output),length(output(1).marginal_likelihood.Rsp.likelihoods));
rsp_value = output(1).marginal_likelihood.Rsp.values;
for i = 1:length(output)
    rsp_lik(i,:) = output(i).marginal_likelihood.Rsp.likelihoods;
end

rsp_lik_heatmap = zeros(length(rsp_value),length(true_di));
for i = 1:cell_num
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_rsp_lik = mean(rsp_lik(ind_lower:ind_upper,:),1);
rsp_lik_heatmap(:,i) = avg_rsp_lik;
end
% extend y axie: should be the number of x (simulate 'true DI' bin)
rsp_lik_pcolor = [rsp_lik_heatmap;zeros(1,length(true_di))];
%extend x axie: should be the number of y (likelihood bin)
rsp_lik_pcolor = [rsp_lik_pcolor,zeros(1+length(rsp_value),1)];
figure(5),
% tiledlayout flow
% nexttile
heatmap(true_di,flip(output(1).marginal_likelihood.Rsp.values),flip(rsp_lik_heatmap),'Colormap',parula)
xlabel('Value of ^\primetrue DI^\prime')
ylabel('Predict of Rsp Value')
title('Distribution of Rsp Probability')

% nexttile
% pcolor([true_di 1],[0 rsp_value],rsp_lik_pcolor)
% xlabel('Value of ^\primetrue DI^\prime')
% ylabel('Predict of Rsp Value')
% title('Distribution of Rsp Probability')
%%
% extract sigma likelihood data in output structure. Transform the data form into
% plotting suitable form
sig_lik = zeros(length(output),length(output(1).marginal_likelihood.sigma.likelihoods));
sig_value = output(1).marginal_likelihood.sigma.values;
for i = 1:length(output)
    sig_lik(i,:) = output(i).marginal_likelihood.sigma.likelihoods;
end

sig_lik_heatmap = zeros(length(sig_value),length(true_di));
for i = 1:cell_num
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_sig_lik = mean(sig_lik(ind_lower:ind_upper,:),1);
sig_lik_heatmap(:,i) = avg_sig_lik;
end
% extend y axie: should be the number of x (simulate 'true DI' bin)
sig_lik_pcolor = [sig_lik_heatmap;zeros(1,length(true_di))];
%extend x axie: should be the number of y (likelihood bin)
sig_lik_pcolor = [sig_lik_pcolor,zeros(1+length(sig_value),1)];
figure(6),
% tiledlayout flow
% nexttile
heatmap(true_di,flip(output(1).marginal_likelihood.sigma.values),flip(sig_lik_heatmap),'Colormap',parula)
xlabel('Value of ^\primetrue DI^\prime')
ylabel('Predict of sig Value')
title('Distribution of sig Probability')

% nexttile
% pcolor([true_di 1],[0 sig_value],sig_lik_pcolor)
% xlabel('Value of ^\primetrue DI^\prime')
% ylabel('Predict of sig Value')
% title('Distribution of sig Probability')
%%
% extract thetapref likelihood data in output structure. Transform the data form into
% plotting suitable form
thetapref_lik = zeros(length(output),length(output(1).marginal_likelihood.theta_pref.likelihoods));
thetapref_value = output(1).marginal_likelihood.theta_pref.values;
for i = 1:length(output)
    thetapref_lik(i,:) = output(i).marginal_likelihood.theta_pref.likelihoods;
end

thetapref_lik_heatmap = zeros(length(thetapref_value),length(true_di));
for i = 1:cell_num
    ind_lower = experiment_num*(i-1)+1;
    ind_upper = experiment_num*i;
avg_thetapref_lik = mean(thetapref_lik(ind_lower:ind_upper,:),1);
thetapref_lik_heatmap(:,i) = avg_thetapref_lik;
end
% extend y axie: should be the number of x (simulate 'true DI' bin)
thetapref_lik_pcolor = [thetapref_lik_heatmap;zeros(1,length(true_di))];
%extend x axie: should be the number of y (likelihood bin)
thetapref_lik_pcolor = [thetapref_lik_pcolor,zeros(1+length(thetapref_value),1)];
figure(7),
% tiledlayout flow
% nexttile
heatmap(true_di,flip(output(1).marginal_likelihood.theta_pref.values),flip(thetapref_lik_heatmap),'Colormap',parula)
xlabel('Value of ^\primetrue DI^\prime')
ylabel('Predict of thetapref Value')
title('Distribution of thetapref Probability')

% nexttile
% pcolor([true_di 1],[0 thetapref_value],thetapref_lik_pcolor)
% xlabel('Value of ^\primetrue DI^\prime')
% ylabel('Predict of thetapref Value')
% title('Distribution of thetapref Probability')
%%
figure(8),
for i = 1:5
subplot(1,5,i),hold on;
plot(data_ideal{i}.angle,data_ideal{i}.responses,'k');
end
