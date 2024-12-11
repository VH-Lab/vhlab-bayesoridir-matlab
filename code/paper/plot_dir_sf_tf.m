function plot_dir_sf_tf(dir_doc)


[dirinfo] = extract_direction_data(dir_doc);

dx = mean(diff(dirinfo(1).angle));
X_total = dirinfo(1).angle(end)-dirinfo(1).angle(1);
dy = 0.1;

% step 1: find the maximum value of the data

M = -Inf;

for i=1:numel(dirinfo),
    M = max(M,max(dirinfo(i).mean_responses));
end;

sfs = unique([dirinfo.sf]);
tfs = unique([dirinfo.tf]);

for i=1:numel(dirinfo),
    s = find(sfs==dirinfo(i).sf);
    t = find(tfs==dirinfo(i).tf);
    hold on;
    plot( (t-1)*(dx+X_total)+dirinfo(i).angle, (s-1)*(dy+1)+dirinfo(i).mean_responses/M );
    ylabel('spacial')
    xlabel('temporal')
end;

