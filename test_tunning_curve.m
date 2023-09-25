ang = 0:359;
size(ang)
rp = 7;
op = 90;
alpha = 0.5;
sig = 180;
rsp = 1;
% one gaussian distribution in angle range [0,360)
vrsp1 = alpha .* rp .* exp(-0.5*angdiff(ang-op+180).^2./sig.^2);
vrsp2 = rp .* exp(-0.5*angdiff(ang-op).^2./sig.^2);
% two gaussian distribution in angle range [0,360)
vrsp3 = rp .* exp(-0.5*angdiff(ang-op).^2./sig.^2) + alpha .* rp .* exp(-0.5*angdiff(ang-op+180).^2./sig.^2);
% two gaussian distribution with a background offset in angle range [0,360)
vrsp4 = rsp + rp .* exp(-0.5*angdiff(ang-op).^2./sig.^2) + alpha .* rp .* exp(-0.5*angdiff(ang-op+180).^2./sig.^2);


tiledlayout(2,2);

% plot three different curves
nexttile,
plot(ang,vrsp1);
ylim([0,8]);
title('prefer direction gaussian distribution');

nexttile,
plot(ang,vrsp2);
ylim([0,8]);
title('null direction gaussian distribution')

nexttile,
plot(ang,vrsp3);
ylim([3,11]);
title('two gaussian distribution');

nexttile,
plot(ang,vrsp4);
ylim([4,12]);
title('two gaussian distribution with background offset');