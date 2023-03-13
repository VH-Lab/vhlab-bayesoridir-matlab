function out = bayes_grid_function(data)
out.data_mean = mean(data,1);
% build 71 x 180 x 71 x 61 grids
out.RP = 85:155;
out.OP = 0:359;
out.RN = 25:95;
out.SIG = 40:100;
out.lik = zeros(length(out.RP),length(out.OP),length(out.RN),length(out.SIG));
ang = 0:359;
%S = median(std(ang,1));
for i = 1:length(out.RP)
    for j = 1:length(out.OP)
        for k = 1:length(out.RN)  
            for h = 1:length(out.SIG)   
                vrsp = out.RP(i) * exp(-0.5*(ang-out.OP(j)).^2/out.SIG(h)^2) + out.RN(k) * exp(-0.5*(ang-(out.OP(j)+180)).^2/out.SIG(h)^2);
                out.lik(i,j,k,h) = exp(-0.5*sum((out.data_mean - vrsp).^2)/5.^2);
            end
        end 
    end
end
