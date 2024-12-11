function hist_index = debug_discretize(debug, debug3)

oi = debug.OI;
di = debug.DI;
cv = debug.CV;
dcv = debug.DCV;

oi_bins  = 0:0.05:1;
di_bins  = 0:0.05:1;
cv_bins  = 0:0.05:1;
dcv_bins = 0:0.05:1;

oi_index = discretize(oi,oi_bins);
di_index = discretize(di,di_bins);
cv_index = discretize(cv,cv_bins);
dcv_index = discretize(dcv,dcv_bins);

oi_index(isnan(oi_index)) = numel(oi_bins)+1;

hist_index = sub2ind(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],oi_index,di_index,cv_index,dcv_index);

IBIG = find(hist_index-double(debug3.histogram_index)==-21296);
 % these all have DI = 0, DCV = 1

ISMALL = find(hist_index-double(debug3.histogram_index)==-10648);
 % these are all a fixed value of DI, DCV, OI, almost fixed CV, only 17 values

[oi_index_alt,di_index_alt,cv_index_alt,dcv_index_alt] = ind2sub(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],...
	hist_index);

[oi_index_alt2,di_index_alt2,cv_index_alt2,dcv_index_alt2] = ind2sub(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],...
	debug3.histogram_index);

oi_index_alt3 = vlt.data.ind2subD(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],debug3.histogram_index,1);
di_index_alt3 = vlt.data.ind2subD(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],debug3.histogram_index,2);
cv_index_alt3 = vlt.data.ind2subD(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],debug3.histogram_index,3);
dcv_index_alt3 = vlt.data.ind2subD(1+[numel(oi_bins) numel(di_bins) numel(cv_bins) numel(dcv_bins)],debug3.histogram_index,4);
