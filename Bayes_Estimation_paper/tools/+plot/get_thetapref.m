function get_thetapref(true_index_value,filename,cell_num,isoi)
% extract thetapref likelihood data in output structure, making a plotting suitable form
%   TRUE_INDEX_VALUE- Different value of underlying OI or DI
%   FILENAME        - Extract bayes likelihood values in output structure
%   CELL_NUM        - The number of different underlying OI or DI
%   ISOI    - true  - true values are from OI
%           - false - true values are from DI
titlename = strrep(filename,'_',' ');
thetapref_lik_heatmap = [];
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
    thetapref_lik = zeros(length(output),length(output(1).marginal_likelihood.theta_pref.likelihoods));
    thetapref_value = output(1).marginal_likelihood.theta_pref.values;
    for ind = 1:length(output)
        thetapref_lik(ind,:) = output(ind).marginal_likelihood.theta_pref.likelihoods;
    end
    heatmap_grid = zeros(length(thetapref_value),length(true_index_value));
    for ind = 1:cell_num
        ind_lower = experiment_num*(ind-1)+1;
        ind_upper = experiment_num*ind;
        avg_thetapref_lik = mean(thetapref_lik(ind_lower:ind_upper,:),1);
        heatmap_grid(:,ind) = avg_thetapref_lik;
    end
    thetapref_lik_heatmap(:,:,i) = heatmap_grid;
end
maxlim = max(thetapref_lik_heatmap,[],'all')
figure(),
tiledlayout flow
for i = 1:length(filename)
    nexttile,
    f1 = heatmap(true_index_value,round(flip(output(1).marginal_likelihood.theta_pref.values),2),flip(thetapref_lik_heatmap(:,:,i)), ...
        'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,maxlim]);
    if isoi
        xlabel('Value of ^\primetrue OI^\prime')
    else
        xlabel('Value of ^\primetrue DI^\prime')
    end
    ylabel('Predict of thetapref Value')
    title([titlename(i) 'Distribution of thetapref Probability'])
end