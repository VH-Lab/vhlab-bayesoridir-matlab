function get_sigma(true_index_value,filename,cell_num,isoi)
% extract sigma likelihood data in output structure, making a plotting suitable form
%   TRUE_INDEX_VALUE- Different value of underlying OI or DI
%   FILENAME        - Extract bayes likelihood values in output structure
%   CELL_NUM        - The number of different underlying OI or DI
%   ISOI    - true  - true values are from OI
%           - false - true values are from DI
titlename = strrep(filename,'_',' ');
sig_lik_heatmap = [];
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
    sig_lik = zeros(length(output),length(output(1).marginal_likelihood.sigma.likelihoods));
    sig_value = output(1).marginal_likelihood.sigma.values;
    for ind = 1:length(output)
        sig_lik(ind,:) = output(ind).marginal_likelihood.sigma.likelihoods;
    end
    heatmap_grid = zeros(length(sig_value),length(true_index_value));
    for ind = 1:cell_num
        ind_lower = experiment_num*(ind-1)+1;
        ind_upper = experiment_num*ind;
        avg_sig_lik = mean(sig_lik(ind_lower:ind_upper,:),1);
        heatmap_grid(:,ind) = avg_sig_lik;
    end
    sig_lik_heatmap(:,:,i) = heatmap_grid;
end
maxlim = max(sig_lik_heatmap,[],'all')
figure(),
tiledlayout(3,4)
for i = 1:length(filename)
    nexttile,
    heatmap(true_index_value,round(flip(output(1).marginal_likelihood.sigma.values),2),flip(sig_lik_heatmap(:,:,i)), ...
        'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,maxlim])
    if isoi
        xlabel('Value of ^\primetrue OI^\prime')
    else
        xlabel('Value of ^\primetrue DI^\prime')
    end
    ylabel('Predict of Sigma Value')
    title([titlename(i) 'Distribution of Sigma Probability'])
end