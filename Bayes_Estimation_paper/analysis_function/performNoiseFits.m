
load sst.mat
load PV.mat
load VIP.mat
load Ex.mat

sst_m = []; sst_v = []; sst_n = 0;

for i=1:numel(sst_meanResponse)
    if ~isempty(sst_meanResponse{i})
        sst_m = sst_meanResponse{i}(:); 
        sst_v = sst_varResponse{i}(:);
        sst_n = sst_n + size(sst_meanResponse{i},2);
    end
end

pv_m = []; pv_v = []; pv_n = 0;

for i=1:numel(PV_meanResponse)
    if ~isempty(PV_meanResponse{i})
        pv_m = PV_meanResponse{i}(:); 
        pv_v = PV_varResponse{i}(:);
        pv_n = pv_n + size(PV_meanResponse{i},2);
    end
end

vip_m = []; vip_v = []; vip_n = 0;

for i=1:numel(VIP_meanResponse)
    if ~isempty(VIP_meanResponse{i})
        vip_m = VIP_meanResponse{i}(:); 
        vip_v = VIP_varResponse{i}(:);
        vip_n = vip_n + size(VIP_meanResponse{i},2);
    end
end

ex_m = []; ex_v = []; ex_n = 0;

for i=1:numel(Ex_meanResponse)
    if ~isempty(Ex_meanResponse{i})
        ex_m = Ex_meanResponse{i}(:); 
        ex_v = Ex_varResponse{i}(:);
        ex_n = ex_n + size(Ex_meanResponse{i},2);
    end
end

subplot(3,2,3);
[sst_noisemodel] = vis.bayes.noise.fit_proportional_noise_plus_c(abs(sst_m),sqrt(sst_v),1);
title(['SST Interneurons (' int2str(sst_n) '): ' mat2str(sst_noisemodel,3)]);
axis([1e-7 1.5 0.001 1]);
axis square
set(gca,'xtick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])

subplot(3,2,4);
[pv_noisemodel] = vis.bayes.noise.fit_proportional_noise_plus_c(abs(pv_m),sqrt(pv_v),1);
title(['PV Interneurons (' int2str(pv_n) '): ' mat2str(pv_noisemodel,3)]);
axis([1e-7 1.5 0.001 1]);
axis square
set(gca,'xtick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])

subplot(3,2,5);
[vip_noisemodel] = vis.bayes.noise.fit_proportional_noise_plus_c(abs(vip_m),sqrt(vip_v),1);
title(['VIP Interneurons (' int2str(vip_n) '): ' mat2str(vip_noisemodel,3)]);
axis([1e-7 1.5 0.001 1]);
axis square
set(gca,'xtick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])

subplot(3,2,6);
[ex_noisemodel] = vis.bayes.noise.fit_proportional_noise_plus_c(abs(ex_m),sqrt(ex_v),1);
title(['EX cells (' int2str(ex_n) '):' mat2str(ex_noisemodel,3)]);
axis([1e-7 1.5 0.001 1]);
axis square
set(gca,'xtick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1])