function get_rpref(true_index_value,filename,cell_num,isoi)
% extract rpref likelihood data in output structure, making a plotting suitable form
%   TRUE_INDEX_VALUE- Different value of underlying OI or DI
%   FILENAME        - Extract bayes likelihood values in output structure
%   CELL_NUM        - The number of different underlying OI or DI
%   ISOI    - true  - true values are from OI
%           - false - true values are from DI
titlename = strrep(filename,'_',' ');
rp_lik_heatmap = [];
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
    rp_lik = zeros(length(output),length(output(1).marginal_likelihood.Rp.likelihoods));
    rp_value = output(1).marginal_likelihood.Rp.values;
    for ind = 1:length(output)
        rp_lik(ind,:) = output(ind).marginal_likelihood.Rp.likelihoods;
    end
    heatmap_grid = zeros(length(rp_value),length(true_index_value));
    for ind = 1:cell_num
        ind_lower = experiment_num*(ind-1)+1;
        ind_upper = experiment_num*ind;
        avg_rp_lik = mean(rp_lik(ind_lower:ind_upper,:),1);
        heatmap_grid(:,ind) = avg_rp_lik;
    end
    rp_lik_heatmap(:,:,i) = heatmap_grid;
end
maxlim = max(rp_lik_heatmap,[],'all')
figure(),
tiledlayout flow
for i = 1:length(filename)
    nexttile,
    heatmap(true_index_value,round(flip(output(1).marginal_likelihood.Rp.values),2),flip(rp_lik_heatmap(:,:,i)), ...
        'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,maxlim])
    if isoi
        xlabel('Value of ^\primetrue OI^\prime')
    else
        xlabel('Value of ^\primetrue DI^\prime')
    end
    ylabel('Predict of Rpref Value')
    title([titlename(i) 'Distribution of Rpref Probability'])
end