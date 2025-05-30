function inPercentile = isInCentralPercentile(values, probabilities, testValue, testPercentiles, options)
%isInCentralPercentile Determines if a value is within the central N% of a distribution.
%
%   inPercentile = isInCentralPercentile(values, probabilities, testValue, testPercentiles)
%   Determines if testValue lies within the central N percent of a probability
%   distribution defined by 'values' and 'probabilities'. The "center" of the
%   distribution is defined by its mode. The distribution is modeled as a
%   piece-wise constant PDF based on the input points.
%
%   inPercentile = isInCentralPercentile(..., Name, Value) specifies options
%   using one or more Name, Value pair arguments.
%
%   Inputs:
%       values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
%           A vector containing the discrete values the random variable can take.
%           Must be sorted in ascending order, or will be sorted internally.
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSize(values, probabilities)}
%           A vector containing the probabilities corresponding to each element
%           in 'values'. The sum of probabilities does not necessarily need to be 1;
%           it will be normalized internally. Must have the same size as 'values'.
%
%       testValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
%           The scalar value to test against the central percentile intervals.
%
%       testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                             mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
%           A vector of percentages (0 to 100). For each percentage N in
%           this vector, the function determines if testValue is within the
%           central N% interval of the distribution.
%
%   Optional Name-Value Pairs:
%       PlotType ('none' | 'pdf' | 'cdf')
%           Specifies whether to generate a plot.
%           'none' (default): No plot is generated.
%           'pdf': Plots the modeled Probability Density Function (PDF), the mode,
%                  the boundaries for each central percentile interval, and the testValue.
%           'cdf': Plots the modeled Cumulative Distribution Function (CDF), the mode,
%                  the boundaries for each central percentile interval, and the testValue.
%
%   Outputs:
%       inPercentile (:,1) logical
%           A logical vector of the same size as 'testPercentiles'. Each element
%           is true if testValue falls within the central interval defined by
%           the corresponding percentage in testPercentiles, and false otherwise.
%
%   Distribution Modeling:
%   The function models the probability distribution as a piece-wise constant
%   Probability Density Function (PDF). Each pair (values(i), probabilities(i))
%   defines the probability mass within a bin. The bin associated with values(i)
%   extends halfway to the previous value point (values(i-1)) and halfway to
%   the next value point (values(i+1)).
%   - Bin edges are calculated at the midpoints between consecutive values.
%   - For the first value, the bin extends from values(1) - (values(2)-values(1))/2
%     to (values(1)+values(2))/2.
%   - For the last value, the bin extends from (values(end-1)+values(end))/2
%     to values(end) + (values(end)-values(end-1))/2.
%   The total probability mass (sum of input 'probabilities') is normalized to 1
%   to form a valid PDF.
%
%   Central Interval Definition:
%   The "center" of the distribution is defined as the mode (the value with the
%   highest associated probability mass). The central N% interval [L, U] is
%   determined such that:
%   1. The mode lies within the interval [L, U].
%   2. The integral of the modeled PDF from L to U equals N/100.
%   3. The interval is defined symmetrically in terms of probability mass
%      around the mode's position within the cumulative distribution. Specifically,
%      if C(x) is the Cumulative Distribution Function (CDF) and m is the mode value,
%      the interval [L, U] satisfies C(U) - C(L) = N/100, and is found by targeting
%      C(L) = C(m) - (N/100)/2 and C(U) = C(m) + (N/100)/2 (clamped to [0, 1]),
%      then finding L and U using the inverse CDF.
%
%   Example:
%       values = [10; 20; 30; 40; 50];
%       probabilities = [0.1; 0.2; 0.4; 0.2; 0.1]; % Mode is 30
%       testVal = 35;
%       percentiles = [50; 90; 95];
%       isInside = isInCentralPercentile(values, probabilities, testVal, percentiles, 'PlotType', 'pdf');
%       % isInside will be [true; true; false] (approx.) and a plot will be generated showing PDF, mode, bounds, and test value.

    arguments
        values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSize(values, probabilities)}
        testValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
        testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                             mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
        options.PlotType {mustBeMember(options.PlotType, {'none', 'pdf', 'cdf'})} = 'none'
    end

    nVal = numel(values);
    nPerc = numel(testPercentiles);
    inPercentile = false(nPerc, 1); % Initialize output
    lower_bounds = NaN(nPerc, 1);   % Initialize storage for bounds
    upper_bounds = NaN(nPerc, 1);   % Initialize storage for bounds


    % --- Input Handling and Preprocessing ---

    % Ensure values are sorted and handle duplicates by summing probabilities
    if ~issorted(values, 'ascend')
        [values, sortIdx] = sort(values);
        probabilities = probabilities(sortIdx);
    end

    [unique_values, ~, ic] = unique(values);
    if numel(unique_values) < nVal
        accumulated_probabilities = accumarray(ic, probabilities, [], @sum);
        values = unique_values;
        probabilities = accumulated_probabilities;
        nVal = numel(values); % Update nVal
        fprintf('Warning: Duplicate values detected. Probabilities were summed.\n');
    end

    % Handle trivial case: single point distribution
    if nVal == 1
        modeValue = values(1); % Mode is the only value
        totalProb = sum(probabilities); % Check total prob even for single point
         if totalProb <= 0
             warning('Total probability mass is zero or negative.');
             inPercentile = false(nPerc, 1);
             return;
         end

        for i = 1:nPerc
           if testPercentiles(i) > 0
               % Any non-zero percentile interval for a single point just contains the point
               lower_bounds(i) = modeValue;
               upper_bounds(i) = modeValue;
               inPercentile(i) = (testValue == modeValue);
           else
               % Zero percentile interval is empty
               lower_bounds(i) = NaN;
               upper_bounds(i) = NaN;
               inPercentile(i) = false;
           end
        end
        % --- Optional Plotting for Single Point ---
        if ~strcmpi(options.PlotType, 'none')
            warning('Plotting is not well-defined for a single-point distribution. Skipping plot.');
        end
        return; % Finished processing single point case
    end

    % Normalize probabilities to sum to 1
    totalProb = sum(probabilities);
    if totalProb <= 0
        warning('Total probability mass is zero or negative.');
        inPercentile = false(nPerc, 1);
        return;
    elseif abs(totalProb - 1.0) > sqrt(eps) % Check if normalization is needed
        normalized_probabilities = probabilities / totalProb;
    else
        normalized_probabilities = probabilities;
    end

    % --- Define PDF Bins ---
    midpoints = (values(1:end-1) + values(2:end)) / 2;
    first_half_width = (midpoints(1) - values(1));
    last_half_width = (values(end) - midpoints(end));
    first_edge = values(1) - max(0, first_half_width);
    last_edge = values(end) + max(0, last_half_width);
    if nVal > 1 && first_edge >= midpoints(1)
        first_edge = values(1) - eps(values(1));
    end
     if nVal > 1 && last_edge <= midpoints(end)
         last_edge = values(end) + eps(values(end));
     end
    edges = [first_edge; midpoints; last_edge];
    widths = diff(edges);
    if any(widths <= 0)
        bad_width_indices = find(widths <= 0);
        problematic_value_indices = unique([bad_width_indices; bad_width_indices+1]);
        problematic_value_indices(problematic_value_indices > nVal) = [];
        warning('Bin calculation resulted in zero or negative widths near values: %s. Check input values.', num2str(values(problematic_value_indices)'));
        inPercentile = false(nPerc, 1);
        return;
    end

    % --- Calculate PDF and CDF ---
    pdf_heights = normalized_probabilities ./ widths;
    cdf_vals = zeros(nVal + 1, 1);
    cdf_vals(2:end) = cumsum(normalized_probabilities);
    cdf_vals(end) = 1.0;

    % --- Find Mode ---
    [~, mode_idx] = max(normalized_probabilities);
    mode_idx = find(normalized_probabilities == normalized_probabilities(mode_idx), 1, 'first');
    modeValue = values(mode_idx); % This is the value corresponding to the mode probability

    % --- Calculate CDF at the Mode Location ---
    cdf_at_bin_start = cdf_vals(mode_idx);
    cdf_at_mode = cdf_at_bin_start + normalized_probabilities(mode_idx) / 2;
    cdf_at_mode = max(0, min(1, cdf_at_mode));

    % --- Calculate Intervals and Test Value ---
    for i = 1:nPerc
        percentile = testPercentiles(i);
        target_mass = percentile / 100.0;

        if target_mass <= 0
            lower_bounds(i) = NaN;
            upper_bounds(i) = NaN;
            inPercentile(i) = false;
            continue;
        elseif target_mass >= 1.0
             lower_bounds(i) = edges(1);
             upper_bounds(i) = edges(end);
        else
            target_cdf_lower = cdf_at_mode - target_mass / 2.0;
            target_cdf_upper = cdf_at_mode + target_mass / 2.0;
            target_cdf_lower = max(0.0, target_cdf_lower);
            target_cdf_upper = min(1.0, target_cdf_upper);

            lower_bounds(i) = interpInverseCDF(target_cdf_lower, edges, cdf_vals, normalized_probabilities, widths, pdf_heights);
            upper_bounds(i) = interpInverseCDF(target_cdf_upper, edges, cdf_vals, normalized_probabilities, widths, pdf_heights);
        end

        % Check if testValue is within the calculated interval
        if ~isnan(lower_bounds(i)) && ~isnan(upper_bounds(i))
             inPercentile(i) = (testValue >= lower_bounds(i)) && (testValue <= upper_bounds(i));
        else
            inPercentile(i) = false;
        end
    end

    % --- Optional Plotting ---
    if ~strcmpi(options.PlotType, 'none')
        figure;
        hold on;
        plotHandles = []; % Store handles for legend
        plotLabels = {};  % Store labels for legend

        if strcmpi(options.PlotType, 'pdf')
            % Plot PDF
            h_pdf = stairs(edges, [pdf_heights; pdf_heights(end)], 'b-', 'LineWidth', 1.5);
            ylabel('Probability Density');
            plotTitle = 'Probability Density Function (PDF)';
            maxY = max(pdf_heights) * 1.1;
             if maxY == 0, maxY = 1; end
             plotHandles(end+1) = h_pdf;
             plotLabels{end+1} = 'PDF';

        elseif strcmpi(options.PlotType, 'cdf')
            % Plot CDF
            h_cdf = plot(edges, cdf_vals, 'b-', 'LineWidth', 1.5);
            ylabel('Cumulative Probability');
            plotTitle = 'Cumulative Distribution Function (CDF)';
            maxY = 1.05;
            plotHandles(end+1) = h_cdf;
            plotLabels{end+1} = 'CDF';
        end

        % Add lines for percentile bounds
        unique_bounds = unique([lower_bounds; upper_bounds]);
        unique_bounds = unique_bounds(~isnan(unique_bounds)); % Remove NaNs

        if ~isempty(unique_bounds)
             bound_style = '--'; % Dashed line for bounds
             bound_color = [0.6 0.6 0.6]; % Grey color
             % Plot all bounds faintly first without labels for clarity
             xline(unique_bounds, bound_style, 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility', 'off');

             % Add one representative line for the legend
             h_bound_leg = xline(NaN, bound_style, 'Color', bound_color, 'LineWidth', 1); % Invisible line for legend entry
             plotHandles(end+1) = h_bound_leg;
             plotLabels{end+1} = 'Central Percentile Bounds';

             % Optionally add labels to bounds (can get crowded)
             % label_strs = arrayfun(@(b) sprintf('Bound (%0.2f)', b), unique_bounds, 'UniformOutput', false);
             % text(unique_bounds + 0.01*(edges(end)-edges(1)), repmat(maxY*0.9, size(unique_bounds)), label_strs, 'Rotation', 90);
        end

        % Add line for the Mode
        h_mode = xline(modeValue, 'k:', 'LineWidth', 1.5); % Black dotted line for mode
        plotHandles(end+1) = h_mode;
        plotLabels{end+1} = sprintf('Mode (%g)', modeValue);

        % Add line for testValue
        h_test = xline(testValue, 'r-.', 'LineWidth', 2); % Red dash-dot line for test value
        plotHandles(end+1) = h_test;
        plotLabels{end+1} = sprintf('Test Value (%g)', testValue);


        xlabel('Value');
        title(plotTitle);
        grid on;
        legend(plotHandles, plotLabels, 'Location', 'best'); % Create legend with specified handles/labels
        % Adjust xlim to show range nicely
        x_range = edges(end) - edges(1);
         if x_range == 0 % Handle case where range is zero
             x_range = abs(edges(1)) * 0.2 + 1; % Add some padding based on value or default
         end
        xlim_low = edges(1) - 0.05 * x_range;
        xlim_high = edges(end) + 0.05 * x_range;
        % Ensure mode and test value are within limits if they fall outside initial range
        xlim_low = min(xlim_low, min(testValue, modeValue) - 0.05 * x_range);
        xlim_high = max(xlim_high, max(testValue, modeValue) + 0.05 * x_range);
        xlim([xlim_low xlim_high]);

        if strcmpi(options.PlotType, 'pdf')
             ylim([0 maxY]);
        else % CDF
             ylim([-0.05 1.05]); % Ensure range 0-1 is visible for CDF
        end

        hold off;
    end

end

% --------------------------------------------------------------------------
% Local Helper Functions
% --------------------------------------------------------------------------

% Custom validation function
function mustHaveSameSize(a, b)
%MUSTHAVESAMESIZE Validate that inputs have the same size.
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSize';
        msg = sprintf('Inputs must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end

% Helper function for Inverse CDF calculation
function x = interpInverseCDF(target_cdf, edges, cdf_vals, normalized_probabilities, widths, pdf_heights)
    % Finds the value x such that CDF(x) = target_cdf.
    nVal = numel(normalized_probabilities);
    tol = sqrt(eps);

    if target_cdf <= 0 + tol, x = edges(1); return; end
    if target_cdf >= 1 - tol, x = edges(end); return; end

    upper_edge_idx = find(cdf_vals >= target_cdf - tol, 1, 'first');

    if abs(cdf_vals(upper_edge_idx) - target_cdf) < tol
        x = edges(upper_edge_idx);
        return;
    end

    bin_idx = upper_edge_idx - 1;

    if bin_idx < 1, bin_idx = 1; end
    if bin_idx > nVal, bin_idx = nVal; end

    cdf_at_bin_start = cdf_vals(bin_idx);
    prob_in_bin = normalized_probabilities(bin_idx);

    if prob_in_bin < tol
        x = edges(bin_idx);
    else
        pdf_h = pdf_heights(bin_idx);
        if pdf_h < tol
            if widths(bin_idx) > tol
                x = edges(bin_idx) + (target_cdf - cdf_at_bin_start) * widths(bin_idx) / prob_in_bin;
            else
                x = edges(bin_idx);
            end
        else
            x = edges(bin_idx) + (target_cdf - cdf_at_bin_start) / pdf_h;
        end
    end
    x = max(edges(bin_idx), min(edges(bin_idx+1), x));
end