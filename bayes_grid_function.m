function [RP,OP,RN,SIG,lik] = bayes_grid_function(stimulus_value,ang)
% out.data_mean = mean(data,1);
% build bayes grid matrix
RP = 0:150;
OP = 0:5:359;
RN = 0:150; 
SIG = 0:50;
Noise_SIG = 0:15;
lik = cell(length(Noise_SIG),1);
for n = 1:length(Noise_SIG)
lik{n} = zeros(length(RP),length(OP),length(RN),length(SIG));
end
% get 5 parameters' posterior probability
for i = 1:length(RP)
    for j = 1:length(OP)
        for k = 1:length(RN)  
            for h = 1:length(SIG)
                vrsp = RP(i) * exp(-0.5*angdiff(ang-OP(j)).^2/SIG(h)^2) + RN(k) * exp(-0.5*angdiff(ang-(OP(j)+180)).^2/SIG(h)^2);
                for n = 1:length(Noise_SIG)
                    prsp = normpdf(vrsp-stimulus_value,0,Noise_SIG(n));
                    lik{n}(i,j,k,h) = prod(prsp);
                end
            end
        end 
    end
end
% TotalLik=[0];
% for n = 1:length(Noise_SIG)
% TotalLik = TotalLik + sum(lik{n},'all');
% end;
% for n = 1:length(Noise_SIG)
%     lik{n} = lik{n}./TotalLik;
% end
%normalized liklihood function

