function get_rsp(true_index_value,filename,cell_num,isoi)
% extract rsp likelihood data in output structure, making a plotting suitable form
%   TRUE_INDEX_VALUE- Different value of underlying OI or DI
%   FILENAME        - Extract bayes likelihood values in output structure
%   CELL_NUM        - The number of different underlying OI or DI
%   ISOI    - true  - true values are from OI
%           - false - true values are from DI
titlename = strrep(filename,'_',' ');
rsp_lik_heatmap = [];
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
    rsp_lik = zeros(length(output),length(output(1).marginal_likelihood.Rsp.likelihoods));
    rsp_value = output(1).marginal_likelihood.Rsp.values;
    for ind = 1:length(output)
        rsp_lik(ind,:) = output(ind).marginal_likelihood.Rsp.likelihoods;
    end
    heatmap_grid = zeros(length(rsp_value),length(true_index_value));
    for ind = 1:cell_num
        ind_lower = experiment_num*(ind-1)+1;
        ind_upper = experiment_num*ind;
        avg_rsp_lik = mean(rsp_lik(ind_lower:ind_upper,:),1);
        heatmap_grid(:,ind) = avg_rsp_lik;
    end
    rsp_lik_heatmap(:,:,i) = heatmap_grid;
end
maxlim = max(rsp_lik_heatmap,[],'all')
figure(),
tiledlayout flow
for i = 1:length(filename)
    nexttile,
    heatmap(true_index_value,round(flip(output(1).marginal_likelihood.Rsp.values),2),flip(rsp_lik_heatmap(:,:,i)), ...
        'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,maxlim])
    if isoi
        xlabel('Value of ^\primetrue OI^\prime')
    else
    xlabel('Value of ^\primetrue DI^\prime')
end
ylabel('Predict of Rsp Value')
title([titlename(i) 'Distribution of Rsp Probability'])
end