function get_alpha(true_index_value,filename,cell_num,isoi)
% extract alpha likelihood data in output structure, making a plotting suitable form
%   TRUE_INDEX_VALUE- Different value of underlying OI or DI
%   FILENAME        - Extract bayes likelihood values in output structure
%   CELL_NUM        - The number of different underlying OI or DI
%   ISOI    - true  - true values are from OI
%           - false - true values are from DI
titlename = strrep(filename,'_',' ');
alpha_lik_heatmap = [];
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
    alpha_lik = zeros(length(output),length(output(1).marginal_likelihood.Alpha.likelihoods));
    alpha_value = output(1).marginal_likelihood.Alpha.values;
    for ind = 1:length(output)
        alpha_lik(ind,:) = output(ind).marginal_likelihood.Alpha.likelihoods;
    end
    heatmap_grid = zeros(length(alpha_value),length(true_index_value));
    for ind = 1:cell_num
        ind_lower = experiment_num*(ind-1)+1;
        ind_upper = experiment_num*ind;
        avg_alpha_lik = mean(alpha_lik(ind_lower:ind_upper,:),1);
        heatmap_grid(:,ind) = avg_alpha_lik;
    end
    alpha_lik_heatmap(:,:,i) = heatmap_grid;
end
maxlim = max(alpha_lik_heatmap,[],'all')
figure()
tiledlayout(3,4)
for i = 1:length(filename)
    nexttile,
    heatmap(true_index_value,round(flip(output(1).marginal_likelihood.Alpha.values),2),flip(alpha_lik_heatmap(:,:,i)), ...
        'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,maxlim]);
    if isoi
        xlabel('Value of ^\primetrue OI^\prime')
    else
        xlabel('Value of ^\primetrue DI^\prime')
    end
    ylabel('Predict of Alpha Value'),
    title([titlename(i) 'Distribution of Alpha Probability']);
end