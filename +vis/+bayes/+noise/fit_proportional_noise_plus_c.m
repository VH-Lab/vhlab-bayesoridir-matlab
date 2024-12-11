function [proportional_noise_plus_c] = fit_proportional_noise_plus_c(m,v, plotit)
% FIT_PROPORTIONAL_NOISE_PLUS_C - fit a noise model STDD = C + K * M^S
%
% [PROPORTIONAL_NOISE] = vis.bayes.noise.fit_proportional_noise_plus_c(M,V, [PLOITIT])
%
% Given mean observations and standard deviations, calculate a linear
% model of M vs. V on a log-log scale (base 10) with an added constant.
%
% Prior papers (Softky Koch 1993 for example) have shown that the standard deviation
% of a response of a neuron V is related to the mean response M on a log-log scale. That is,
%   LOG10(V) ~= SLOPE * LOG10(M) + B.
% If we solve for V as a function of M, then 
%   V(M) ~= 10.^B * M^SLOPE.
%
% Here we let the relation have a constant term so that V(M) isn't 0. It doesn't make sense
% for a cell to be noiseless if its average response is 0.  So we solve for 
% 
% V(M) = C + K * M ^ S
%
% PORTIONAL_NOISE is a 3 element vector [C K S].
%
% If PLOTIT is provided and is 1, the data and the fit are plotted in the current axes.
%
% If any M are less than 0, the absolute value is used.
%

if nargin<3,
	plotit = 0;
end;

myfit = fittype('a+b*x^c');
fo = fitoptions(myfit),

fo.StartPoint = [0 1 1];
fo.Lower = [0 -Inf -Inf ];
fo.Upper = [Inf Inf Inf ];
myfit = setoptions(myfit,fo);

[myanswer,gof] = fit(abs(m(:)),v(:),myfit);

proportional_noise_plus_c = [myanswer.a ; myanswer.b; myanswer.c];

if plotit,
	fit_values_x = logspace(log10(0.0001),log10(max(abs(m))),100);
	fit_values_y = myanswer(fit_values_x);
	plot(abs(m),v,'ko');
	hold on;
	plot(fit_values_x,fit_values_y,'b-','linewidth',2);
	set(gca,'xscale','log','yscale','log');
	ylabel('Response standard deviation');
	xlabel('Response mean');
	box off;
end;


