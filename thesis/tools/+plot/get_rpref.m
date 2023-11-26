function get_rpref(true_oi,filename,cell_num,isoi)
% extract rpref likelihood data in output structure, making a plotting suitable form
%
% ISOI  - true  - true values are from OI
%       - false - true values are from DI

titlename = strrep(filename,'_',' ');
figure(),
tiledlayout(3,4)
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
rp_lik = zeros(length(output),length(output(1).marginal_likelihood.Rp.likelihoods));
rp_value = output(1).marginal_likelihood.Rp.values;
for ind = 1:length(output)
    rp_lik(ind,:) = output(ind).marginal_likelihood.Rp.likelihoods;
end
rp_lik_heatmap = zeros(length(rp_value),length(true_oi));
for ind = 1:cell_num
    ind_lower = experiment_num*(ind-1)+1;
    ind_upper = experiment_num*ind;
avg_rp_lik = mean(rp_lik(ind_lower:ind_upper,:),1);
rp_lik_heatmap(:,ind) = avg_rp_lik;
end
% max(rp_lik_heatmap,[],'all')
nexttile
heatmap(true_oi,round(flip(output(1).marginal_likelihood.Rp.values),2),flip(rp_lik_heatmap), ...
    'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,0.15])
if isoi
    xlabel('Value of ^\primetrue OI^\prime')
else
    xlabel('Value of ^\primetrue DI^\prime')
end
ylabel('Predict of Rpref Value')
title([titlename(i) 'Distribution of Rpref Probability'])
end