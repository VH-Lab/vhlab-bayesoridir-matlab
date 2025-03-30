function [boundary,in_boundary] = in_boundaries(values,likelihoods,true_value,bound_percentage)

in_boundary = zeros(size(bound_percentage));

cdf = cumsum(likelihoods);
[~,idx_mode] = max(likelihoods);
cdf_mode = cdf(idx_mode);


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

idx_lower = in_boundary;
idx_upper = in_boundary;
lower_value = in_boundary;
upper_value = in_boundary;

for i = 1:numel(bound_percentage)
    idx_lower(i) = find(cdf > lower_bound(i),1,"first");
    idx_upper(i) = find(cdf < upper_bound(i),1,"last");
    lower_value(i) = values(idx_lower(i));
    upper_value(i) = values(idx_upper(i));
    if true_value >= lower_value(i) && true_value <= upper_value(i)

        in_boundary(i) = 1;

    end
end

boundary.lower_boundaries = lower_bound;
boundary.upper_boundaries = upper_bound;
boundary.lower_value = lower_value;
boundary.upper_value = upper_value;



end