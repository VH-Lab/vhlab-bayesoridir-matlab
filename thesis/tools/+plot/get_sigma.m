function get_sigma(true_oi,filename,cell_num,isoi)
% extract sigma likelihood data in output structure, making a plotting suitable form
%
titlename = strrep(filename,'_',' ');
figure(),
tiledlayout(3,4)
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
sig_lik = zeros(length(output),length(output(1).marginal_likelihood.sigma.likelihoods));
sig_value = output(1).marginal_likelihood.sigma.values;
for ind = 1:length(output)
    sig_lik(ind,:) = output(ind).marginal_likelihood.sigma.likelihoods;
end
sig_lik_heatmap = zeros(length(sig_value),length(true_oi));
for ind = 1:cell_num
    ind_lower = experiment_num*(ind-1)+1;
    ind_upper = experiment_num*ind;
avg_sig_lik = mean(sig_lik(ind_lower:ind_upper,:),1);
sig_lik_heatmap(:,ind) = avg_sig_lik;
end
% max(sig_lik_heatmap,[],'all')
nexttile
heatmap(true_oi,flip(output(1).marginal_likelihood.sigma.values),flip(sig_lik_heatmap), ...
    'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,0.45])
if isoi
    xlabel('Value of ^\primetrue OI^\prime')
else
    xlabel('Value of ^\primetrue DI^\prime')
end
ylabel('Predict of Sigma Value')
title([titlename(i) 'Distribution of Sigma Probability'])
end