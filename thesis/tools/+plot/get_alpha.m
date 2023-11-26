function get_alpha(true_oi,filename,cell_num,isoi)
%
titlename = strrep(filename,'_',' ');
figure()
tiledlayout(3,4)
for i = 1:length(filename)
    load(filename{i}),
    experiment_num = length(output)/cell_num;
    % extract alpha likelihood data in output structure, making a plotting suitable form
    alpha_lik = zeros(length(output),length(output(1).marginal_likelihood.Alpha.likelihoods));
    alpha_value = output(1).marginal_likelihood.Alpha.values;
    for ind = 1:length(output)
        alpha_lik(ind,:) = output(ind).marginal_likelihood.Alpha.likelihoods;
    end
    alpha_lik_heatmap = zeros(length(alpha_value),length(true_oi));
    for ind = 1:cell_num
        ind_lower = experiment_num*(ind-1)+1;
        ind_upper = experiment_num*ind;
        avg_alpha_lik = mean(alpha_lik(ind_lower:ind_upper,:),1);
        alpha_lik_heatmap(:,ind) = avg_alpha_lik;
    end
    nexttile,
    heatmap(true_oi,round(flip(output(1).marginal_likelihood.Alpha.values),2),flip(alpha_lik_heatmap), ...
        'Colormap',parula,'CellLabelColor','none','ColorLimits',[0,0.25]);
if isoi
    xlabel('Value of ^\primetrue OI^\prime')
else
    xlabel('Value of ^\primetrue DI^\prime')
end
    ylabel('Predict of Alpha Value'),
    title([titlename(i) 'Distribution of Alpha Probability']);
end