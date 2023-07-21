function Probability = bayes_grid_function_fixed_noise(I,stimulus_value,ang,noise_mdl)
% extract coefficients of noise linear regression model
offset = noise_mdl(1);
slope = noise_mdl(2);
% build bayes grid matrix
Probability.Rp = zeros(length(I.Rp));
Probability.Op = zeros(length(I.Op));
Probability.Alpha = zeros(length(I.Alpha));
Probability.Sig = zeros(length(I.Sig));
Probability.Lik = zeros(length(I.Rp),length(I.Op),length(I.Alpha),length(I.Sig));

% get 5 parameters' posterior probability(liklihood)
for i = 1:length(I.Rp)
    i
    for j = 1:length(I.Op)
        for k = 1:length(I.Alpha)
            for h = 1:length(I.Sig)
                vrsp = I.Rp(i) * exp(-0.5*angdiff(ang-I.Op(j)).^2/I.Sig(h)^2) + I.Alpha(k)* I.Rp(i) * exp(-0.5*angdiff(ang-(I.Op(j)+180)).^2/I.Sig(h)^2);
                prsp = normpdf(vrsp-stimulus_value,zeros(size(vrsp)),10.^offset*vrsp.^slope);
                Probability.Lik(i,j,k,h) = prod(prsp);

            end
        end 
    end
end
%normalized liklihood function
Probability.Rp =sum(sum(sum(Probability.Lik,4),3),2);
Probability.Rp = Probability.Rp/sum(Probability.Rp); 
Probability.Op = sum(sum(sum(Probability.Lik,4),3),1)';
Probability.Op = Probability.Op/sum(Probability.Op);
Probability.Alpha = squeeze(sum(sum(sum(Probability.Lik,4),2),1)); 
Probability.Alpha = Probability.Alpha/sum(Probability.Alpha);
Probability.Sig = squeeze(sum(sum(sum(Probability.Lik,3),2),1));
Probability.Sig = Probability.Sig/sum(Probability.Sig);

% TotalLik=[0];
% for n = 1:length(Noise_SIG)
% TotalLik = TotalLik + sum(lik{n},'all');
% end;
% for n = 1:length(Noise_SIG)
%     lik{n} = lik{n}./TotalLik;
% end

% draw maximum lik response function v.s. raw data
[M,ind] = max(Probability.Lik,[],'all');
[i,j,k,h] = ind2sub(size(Probability.Lik),ind);
vrsp = I.Rp(i) * exp(-0.5*angdiff(ang-I.Op(j)).^2/I.Sig(h)^2) + I.Alpha(k)* I.Rp(i) * exp(-0.5*angdiff(ang-(I.Op(j)+180)).^2/I.Sig(h)^2);
plot(ang,stimulus_value,'k',ang,vrsp,'b');
end