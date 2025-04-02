function [boundaries,in_boundary] = in_boundaries(values,likelihoods,true_value,bound_percentage)

middle_edges = (values(1:end-1) + values(2:end)) / 2;
lower_edge = values(1);
upper_edge = values(end);
bin_edges = [lower_edge;middle_edges(:);upper_edge];
bin_widths = diff(bin_edges);

% figure();
% for i = 1:length(likelihoods)
%     x = bin_edges(i);
%     w = bin_widths(i);
%     h = likelihoods(i);
%     rectangle('Position', [x, 0, w, h], 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k');
% end
% hold on,
% plot(values,likelihoods,'k');
% ylim([0, max(likelihoods)*1.1]);
% hold off;

mass = likelihoods .* bin_widths;
cdf = [0;cumsum(mass)];%adding cdf start point
cdf = cdf/cdf(end);
in_boundary = zeros(size(bound_percentage));
[~,idx_mode] = max(likelihoods);
mode = values(idx_mode);
cdf_mode = interp1(bin_edges,cdf,mode,'linear');

% Define offsets for the central intervals
offset = bound_percentage./2;

% Calculate rp interval boundaries
% Include the boundary effect when the lower or upper interval boundary
% hit the edge of cdf.
lower_bound = cdf_mode - offset;
upper_bound = cdf_mode + offset;
lower_bound(cdf_mode - offset < 0) = 0;
upper_bound(cdf_mode - offset < 0) = 2*offset(cdf_mode - offset < 0);
upper_bound(cdf_mode + offset > 1) = 1;
lower_bound(cdf_mode + offset > 1) = 1 - 2*offset(cdf_mode + offset > 1);

[cdf_edge,idx_edge] = unique(cdf);
bin_edges = bin_edges(idx_edge);

lower_value = interp1(cdf_edge,bin_edges,lower_bound,'linear');
upper_value = interp1(cdf_edge,bin_edges,upper_bound,'linear');

for i = 1:numel(bound_percentage)
    if true_value >= lower_value(i) && true_value <= upper_value(i)

        in_boundary(i) = 1;

    end
end

boundaries.lower_boundaries = lower_bound;
boundaries.upper_boundaries = upper_bound;
boundaries.lower_value = lower_value;
boundaries.upper_value = upper_value;

% figure(),
% plot(bin_edges,cdf_edge,'r');

end