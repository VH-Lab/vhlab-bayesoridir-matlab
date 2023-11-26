function get_rsp(true_oi,filename,cell_num,isoi)
% extract rsp likelihood data in output structure, making a plotting suitable form
%
titlename = strrep(filename,'_',' ');
figure(),
tiledlayout(3,4)
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
rsp_lik = zeros(length(output),length(output(1).marginal_likelihood.Rsp.likelihoods));
rsp_value = output(1).marginal_likelihood.Rsp.values;
for ind = 1:length(output)
    rsp_lik(ind,:) = output(ind).marginal_likelihood.Rsp.likelihoods;
end
rsp_lik_heatmap = zeros(length(rsp_value),length(true_oi));
for ind = 1:cell_num
    ind_lower = experiment_num*(ind-1)+1;
    ind_upper = experiment_num*ind;
avg_rsp_lik = mean(rsp_lik(ind_lower:ind_upper,:),1);
rsp_lik_heatmap(:,ind) = avg_rsp_lik;
end
% max(rsp_lik_heatmap,[],'all')
nexttile
heatmap(true_oi,round(flip(output(1).marginal_likelihood.Rsp.values),2),flip(rsp_lik_heatmap), ...
    'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,1])
if isoi
    xlabel('Value of ^\primetrue OI^\prime')
else
    xlabel('Value of ^\primetrue DI^\prime')
end
ylabel('Predict of Rsp Value')
title([titlename(i) 'Distribution of Rsp Probability'])
end