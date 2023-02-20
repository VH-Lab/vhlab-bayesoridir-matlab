 function out = my_function(x)
out.y = sin(x);
out.x = x;
out(2).y = cos(x);
out(2).x = 2*x;
%% 
clc
clear
sig = 4;
x = -15:0.01:15
y = 500*exp(-0.5*x.^2/sig.^2)
f1 = figure,
histogram(sig*randn(1,10000))
hold on
plot(x,y,'r')

