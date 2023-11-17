function [proportional_noise] = fit_proportional_noise(m,v, plotit)
% FIT_PROPORTIONAL_NOISE
%
% [PROPORTIONAL_NOISE] = vis.bayes.noise.fit_proportional_noise(M,V)
%
% Given mean observations and standard deviations, calculate a linear
% model of M vs. V on a log-log scale (base 10).
%
% 

if nargin<3,
    plotit = 0;
end;

index = find(m>0 & v>0);

m2 = m(index);
v2 = v(index);

md1 = fitlm(log10(m2),log10(v2));

proportional_noise = md1.Coefficients{:,1};

if plotit,
  % plot m2, v2, and the fit
    offset = proportional_noise(1);
    slope = proportional_noise(2);
    x = logspace(-3,4,100);
    y = 10.^offset*x.^slope;
    plot(x,y,'o')
    set(gca,'xscale','log','yscale','log')
end;