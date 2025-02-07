clear;close;clc;

% Supplementary Figure 1 — An example of the joint marginal likelihood
% distribution.
load my_fig1.mat

value_rp = bayes_w.marginal_likelihood.Rp.values;
value_alpha = bayes_w.marginal_likelihood.Alpha.values;
value_theta = bayes_w.marginal_likelihood.theta_pref.values;
value_sigma = bayes_w.marginal_likelihood.sigma.values;
value_rsp = bayes_w.marginal_likelihood.Rsp.values;

figure(1);
hold on,
plot(ideal_w.angle,ideal_w.responses,'k','LineWidth',1),% show well tuned ideal curve
errorbar(noisy_w.angle,noisy_w.mean_responses,noisy_w.responses_stderr,'b*'),%show well tuned noisy points
plot(0:359,bayes_w.maximum_likelihood.parameters.tunning_curve,'r','LineWidth',1),
xlim([-5,365]),
ylim([0,20]),
yticks(0:4:20),
set(gca,'FontSize',8),
legend('Ideal','+50% Noise','Bayes Estimation','FontSize',7),
xlabel('Direction of Stimuli Motion','FontSize',10),
ylabel('Response(Hz)','FontSize',10),
title('A','FontSize',12);

figure(2);
h = heatmap(rp_thetapref);
Xvalue = value_theta;
Xnum = numel(Xvalue);
Xlabels = cell(Xnum, 1);
for i = 1:Xnum
    if mod(Xvalue(i), 50) == 0
        Xlabels{i} = num2str(Xvalue(i));
    else
        Xlabels{i} = '';
    end
end
h.XDisplayLabels = Xlabels;

Yvalue = zeros(size(value_rp));
Yvalue([1 10 20 30 40 50 60]) = 1;
Yvalue = Yvalue .* value_rp;
Ynum = numel(Yvalue);
Ylabels = cell(Ynum, 1);
for i = 1:Ynum
    if Yvalue(i) ~= 0
        Ylabels{i} = num2str(Yvalue(i));
    else
        Ylabels{i} = '';
    end
end
h.YDisplayLabels = Ylabels;
h.XLabel = '\theta_{pref}';
h.YLabel = 'R_{p}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';

figure(3);
h = heatmap(rp_alpha);
Xvalue = value_alpha;
Xnum = numel(Xvalue);
Xlabels = cell(Xnum, 1);
for i = 1:Xnum
    if mod(Xvalue(i), 0.1) == 0
        Xlabels{i} = num2str(Xvalue(i));
    else
        Xlabels{i} = '';
    end
end
h.XDisplayLabels = Xlabels;
Yvalue = zeros(size(value_rp));
Yvalue([1 10 20 30 40 50 60]) = 1;
Yvalue = Yvalue .* value_rp;
Ynum = numel(Yvalue);
Ylabels = cell(Ynum, 1);
for i = 1:Ynum
    if Yvalue(i) ~= 0
        Ylabels{i} = num2str(Yvalue(i));
    else
        Ylabels{i} = '';
    end
end
h.YDisplayLabels = Ylabels;
h.XLabel = '\alpha';
h.YLabel = 'R_{p}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';

figure(4);
h = heatmap(rp_sigma);
Xvalue = value_sigma;
Xnum = numel(Xvalue);
Xlabels = cell(Xnum, 1);
for i = 1:Xnum
    if mod(Xvalue(i), 5) == 0
        Xlabels{i} = num2str(Xvalue(i));
    else
        Xlabels{i} = '';
    end
end
h.XDisplayLabels = Xlabels;
Yvalue = zeros(size(value_rp));
Yvalue([1 10 20 30 40 50 60]) = 1;
Yvalue = Yvalue .* value_rp;
Ynum = numel(Yvalue);
Ylabels = cell(Ynum, 1);
for i = 1:Ynum
    if Yvalue(i) ~= 0
        Ylabels{i} = num2str(Yvalue(i));
    else
        Ylabels{i} = '';
    end
end
h.YDisplayLabels = Ylabels;
h.XLabel = '\sigma';
h.YLabel = 'R_{p}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';

figure(5);
h = heatmap(thetapref_rsp);
Xvalue = zeros(size(value_rsp));
Xvalue([1 10 20 30 40 50 60]) = 1;
Xvalue = Xvalue .* value_rsp;
Xnum = numel(Xvalue);
Xlabels = cell(Xnum, 1);
for i = 1:Xnum
    if Xvalue(i) ~= 0
        Xlabels{i} = num2str(Xvalue(i));
    else
        Xlabels{i} = '';
    end
end
h.XDisplayLabels = Xlabels;
Yvalue = value_theta;
Ynum = numel(Yvalue);
Ylabels = cell(Ynum, 1);
for i = 1:Ynum
    if mod(Yvalue(i), 50) == 0
        Ylabels{i} = num2str(Yvalue(i));
    else
        Ylabels{i} = '';
    end
end
h.YDisplayLabels = Ylabels;
h.XLabel = 'Rsp';
h.YLabel = '\theta_{pref}';
h.Colormap = parula;
h.GridVisible = 'off';
h.CellLabelColor = 'none';
%%
% export figures
supF1a = figure(2);
supF1b = figure(3);
supF1c = figure(4);
supF1d = figure(5);


path = 'D:\GitHub\vhlab-bayesoridir-matlab\Bayes_Estimation_paper\figures\noise_mdl\kcs\';

exportgraphics(supF1a,[path 'supplementary_Figure_1a.pdf'],"ContentType","vector")
exportgraphics(supF1b,[path 'supplementary_Figure_1b.pdf'],"ContentType","vector")
exportgraphics(supF1c,[path 'supplementary_Figure_1c.pdf'],"ContentType","vector")
exportgraphics(supF1d,[path 'supplementary_Figure_1d.pdf'],"ContentType","vector")
