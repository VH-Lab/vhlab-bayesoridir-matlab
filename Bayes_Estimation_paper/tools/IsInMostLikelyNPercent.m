function isInside = IsInMostLikelyNPercent(values, probabilities, testValue, testPercentiles, options)
%IsInMostLikelyNPercent Determines if a value is within the most likely N% of a distribution.
%
%   isInside = IsInMostLikelyNPercent(values, probabilities, testValue, testPercentiles)
%   Determines if testValue lies within the region(s) encompassing the most
%   likely N percent of a probability distribution defined by 'values' and
%   'probabilities'. The distribution is modeled as a piece-wise constant PDF.
%
%   Inputs:
%       values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
%           A vector containing the discrete values the random variable can take.
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSize(values, probabilities)}
%           A vector containing the probabilities corresponding to each element
%           in 'values'. Will be normalized internally.
%
%       testValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
%           The scalar value to test.
%
%       testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                             mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
%           A vector of percentages (0 to 100). For each percentage N, the
%           function determines if testValue is within the most likely N% region.
%
%   Optional Name-Value Pairs:
%       Plot ('none' | 'pdf')
%           Specifies whether to generate a plot.
%           'none' (default): No plot is generated.
%           'pdf': Plots the modeled Probability Density Function (PDF),
%                  horizontal lines indicating the probability density threshold
%                  for each specified percentile, and the testValue.
%
%   Outputs:
%       isInside (:,1) logical
%           A logical vector of the same size as 'testPercentiles'. Each element
%           is true if testValue falls within the most likely N% region defined
%           by the corresponding percentage in testPercentiles, and false otherwise.
%
%   Definition of Most Likely N% Region:
%   The function first constructs a piece-wise constant PDF from the input
%   values and probabilities (normalized to integrate to 1). It then identifies
%   the unique PDF height values present. To find the region for a given N%,
%   it determines a threshold probability density level 'h_threshold'. This
%   threshold is the highest density value such that the total probability mass
%   contained in all bins with PDF height greater than or equal to 'h_threshold'
%   is at least N/100. The most likely N% region then consists of all value
%   intervals corresponding to these bins (bins where PDF height >= h_threshold).
%   This region may be discontinuous if the PDF is multimodal. The function
%   then checks if testValue falls within any of these intervals.
%
%   Distribution Modeling:
%   (Same as isInCentralPercentile: piece-wise constant PDF based on bins
%   extending halfway between value points, normalized).
%
%   Example 1: Unimodal distribution
%       values = [10; 20; 30; 40; 50];
%       probs = [0.1; 0.2; 0.4; 0.2; 0.1]; % PDF peak at 30
%       testVal = 35; % Near the peak
%       percentiles = [50; 90]; % Check 50% and 90% most likely regions
%       isInside = IsInMostLikelyNPercent(values, probs, testVal, percentiles, 'Plot', 'pdf');
%       % Expected: For N=50, threshold might only include bin around 30.
%       % For N=90, threshold likely lower, includes bins around 20, 30, 40.
%       % testVal=35 is in bin [30+(40-30)/2]=35 boundary? Bin 40 is [35, 45]. PDF height=0.2/10=0.02.
%       % Bin 30 is [25, 35]. PDF height=0.4/10=0.04. 35 is on edge.
%       % isInside likely [true; true] (assuming edges are inclusive).
%
%   Example 2: Bimodal distribution
%       values = [10; 20; 30; 40; 50; 60];
%       probs = [0.3; 0.1; 0.05; 0.05; 0.1; 0.3]; % Peaks near 10 and 60
%       testVal = 25; % Between peaks
%       percentiles = [50; 70; 95];
%       isInside = IsInMostLikelyNPercent(values, probs, testVal, percentiles);
%       % Expected: For N=50, threshold likely includes only bins near 10 and 60.
%       % testVal=25 is likely outside. As N increases, threshold drops, might include middle bins.
%       % isInside might be [false; false; true] (approx).

    arguments
        values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSize(values, probabilities)}
        testValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
        testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                             mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
        options.Plot {mustBeMember(options.Plot, {'none', 'pdf'})} = 'none'
    end

    nVal = numel(values);
    nPerc = numel(testPercentiles);
    isInside = false(nPerc, 1); % Initialize output
    tol = sqrt(eps); % Tolerance for floating point comparisons

    % --- Input Handling and Preprocessing ---
    if isempty(values)
        warning('Input values are empty.');
        return; % Output is already initialized to false
    end
    if ~issorted(values, 'ascend')
        [values, sortIdx] = sort(values);
        probabilities = probabilities(sortIdx);
    end
    [unique_values, ~, ic] = unique(values);
    if numel(unique_values) < nVal
        accumulated_probabilities = accumarray(ic, probabilities, [], @sum);
        values = unique_values;
        probabilities = accumulated_probabilities;
        nVal = numel(values);
        fprintf('Warning: Duplicate values detected. Probabilities were summed.\n');
    end

    % Handle trivial case: single point distribution
    if nVal == 1
        totalProb = sum(probabilities);
        if totalProb <= tol
             warning('Total probability mass is zero or negative.');
             return;
        end
        % For single point, any N>0 includes the point if prob>0.
        for i = 1:nPerc
            isInside(i) = (testPercentiles(i) > tol) && (abs(testValue - values(1)) < tol);
        end
         if ~strcmpi(options.Plot, 'none')
             warning('Plotting is not applicable for a single-point distribution. Skipping plot.');
         end
        return;
    end

    % Normalize probabilities
    totalProb = sum(probabilities);
    if totalProb <= tol
        warning('Total probability mass is zero or negative.');
        return;
    elseif abs(totalProb - 1.0) > tol
        normalized_probabilities = probabilities / totalProb;
    else
        normalized_probabilities = probabilities;
    end

    % --- Define PDF Bins ---
    midpoints = (values(1:end-1) + values(2:end)) / 2;
    first_half_width = max(0, midpoints(1) - values(1)); % Ensure non-negative
    last_half_width = max(0, values(end) - midpoints(end)); % Ensure non-negative
    first_edge = values(1);
    last_edge = values(end);
    if nVal > 1 && abs(first_edge - midpoints(1)) < tol % Handle very close points
        first_edge = values(1) - eps(values(1));
    end
     if nVal > 1 && abs(last_edge - midpoints(end)) < tol
         last_edge = values(end) + eps(values(end));
     end
    edges = [first_edge; midpoints; last_edge];
    widths = diff(edges);
    if any(widths <= tol)
        bad_width_indices = find(widths <= tol);
        problematic_value_indices = unique([bad_width_indices; bad_width_indices+1]);
        problematic_value_indices(problematic_value_indices > nVal) = [];
        warning('Bin calculation resulted in zero or negative widths near values: %s. Check input values.', num2str(values(problematic_value_indices)'));
        return; % Cannot proceed with zero-width bins
    end

    % --- Calculate PDF Heights ---
    pdf_heights = normalized_probabilities ./ widths;
    % Handle potential division by zero if width check failed or probs are huge
    pdf_heights(~isfinite(pdf_heights)) = 0; % Replace Inf/NaN with 0

    % --- Find Threshold Heights for All Percentiles ---
    unique_heights = unique(pdf_heights);
    % Remove effectively zero heights if they exist, check >= h later handles it
    unique_heights = unique_heights(unique_heights > tol);
    sorted_heights = sort(unique_heights, 'descend');
    n_unique_heights = numel(sorted_heights);
    cumulative_mass_at_height = zeros(n_unique_heights, 1);
    bin_indices_for_height = cell(n_unique_heights, 1);

    if isempty(sorted_heights) % All pdf heights were <= tol
        warning('All calculated PDF heights are zero or negligible.');
        % Thresholds will all be 0, isInside will remain false unless N=0
        h_thresholds = zeros(nPerc, 1); % All thresholds are effectively 0
    else
        for j = 1:n_unique_heights
            h = sorted_heights(j);
            % Find bins with height >= current height level h
            % Use >= h - tol for robustness
            current_indices = find(pdf_heights >= h - tol);
            cumulative_mass_at_height(j) = sum(normalized_probabilities(current_indices));
            bin_indices_for_height{j} = current_indices; % Store indices for later use
        end

        % Determine threshold height for each test percentile
        h_thresholds = zeros(nPerc, 1);
        for p_idx = 1:nPerc
            target_mass = testPercentiles(p_idx) / 100.0;
            if target_mass <= tol
                h_thresholds(p_idx) = inf; % Effectively selects no bins
                continue;
            end
            % Find first height level where cumulative mass meets or exceeds target
            threshold_level_idx = find(cumulative_mass_at_height >= target_mass - tol, 1, 'first');

            if isempty(threshold_level_idx)
                % Target mass > 1 (or close to 1 and couldn't reach) -> include all bins
                h_thresholds(p_idx) = 0; % Threshold includes everything non-zero
            else
                h_thresholds(p_idx) = sorted_heights(threshold_level_idx);
            end
        end
    end


    % --- Check Test Value Against Thresholded Regions ---
    for p_idx = 1:nPerc
        current_h_threshold = h_thresholds(p_idx);

        if isinf(current_h_threshold) % Case N=0
             isInside(p_idx) = false;
             continue;
        end

        % Find all bins with PDF height >= threshold
        threshold_bin_indices = find(pdf_heights >= current_h_threshold - tol);

        if isempty(threshold_bin_indices)
            isInside(p_idx) = false; % No bins meet the threshold
            continue;
        end

        % Check if testValue is within any of the intervals defined by these bins
        is_in_any_interval = false;
        for k = 1:numel(threshold_bin_indices)
            bin_idx = threshold_bin_indices(k);
            % Check interval [edges(bin_idx), edges(bin_idx+1)]
            % Use tolerance for boundary checks
            if (testValue >= edges(bin_idx) - tol) && (testValue <= edges(bin_idx+1) + tol)
                is_in_any_interval = true;
                break; % Found it in one interval, no need to check others
            end
        end
        isInside(p_idx) = is_in_any_interval;
    end

    % --- Optional Plotting ---
    if strcmpi(options.Plot, 'pdf')
        figure;
        hold on;
        plotHandles = [];
        plotLabels = {};

        % Plot PDF
        h_pdf = stairs(edges, [pdf_heights; pdf_heights(end)], 'b-', 'LineWidth', 1.5);
        plotHandles(end+1) = h_pdf;
        plotLabels{end+1} = 'PDF';

        % Plot Threshold Lines
        colors = lines(nPerc); % Distinct colors for thresholds
        unique_thresholds = unique(h_thresholds); % Plot each unique threshold once
        plotted_thresh_legend = false;

        for p_idx = 1:nPerc % Iterate through original percentiles to use correct colors
            h_thresh = h_thresholds(p_idx);
             if isinf(h_thresh), continue; end % Don't plot threshold for N=0

             label = sprintf('N=%g%% Threshold (h=%.2g)', testPercentiles(p_idx), h_thresh);
             h_line = yline(h_thresh, '--', label, 'Color', colors(p_idx,:), 'LineWidth', 1);
             h_line.LabelHorizontalAlignment = 'left';

             % Manage legend entries (show only one representative threshold line perhaps?)
             if ~plotted_thresh_legend && h_thresh > tol % Avoid legend entry for threshold=0 unless it's the only one
                  plotHandles(end+1) = h_line;
                  plotLabels{end+1} = 'N% PDF Thresholds';
                  plotted_thresh_legend = true; % Only add one entry for thresholds
             end
             % Ensure subsequent lines for same threshold don't overwrite legend
             h_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
         % If we never added a threshold legend entry (e.g., all thresholds were 0)
         if ~plotted_thresh_legend && any(abs(h_thresholds) < tol & testPercentiles > tol)
              h_line = yline(0, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
              plotHandles(end+1) = h_line;
              plotLabels{end+1} = 'N% PDF Thresholds (at 0)';
              h_line.Annotation.LegendInformation.IconDisplayStyle = 'off';
         end


        % Plot Test Value
        h_test = xline(testValue, 'r-.', 'LineWidth', 2);
        plotHandles(end+1) = h_test;
        plotLabels{end+1} = sprintf('Test Value (%.2g)', testValue);


        ylabel('Probability Density');
        xlabel('Value');
        title('PDF with Most Likely N% Thresholds');
        grid on;
        legend(plotHandles, plotLabels, 'Location', 'best');
        maxY = max(pdf_heights) * 1.1;
        if maxY < tol, maxY = 1; end % Handle case of all zero heights
        ylim([0 maxY]);
         x_range = edges(end) - edges(1);
         if abs(x_range) < tol, x_range = abs(values(1))*0.2+1; end % Handle zero range
         xlim_low = edges(1) - 0.05 * x_range;
         xlim_high = edges(end) + 0.05 * x_range;
         xlim_low = min(xlim_low, testValue - 0.05*x_range); % Ensure test value visible
         xlim_high = max(xlim_high, testValue + 0.05*x_range);% Ensure test value visible
         xlim([xlim_low xlim_high]);
        hold off;
    end

end

% --------------------------------------------------------------------------
% Local Helper / Validation Functions (if not on path)
% --------------------------------------------------------------------------
function mustHaveSameSize(a, b)
%MUSTHAVESAMESIZE Validate that inputs have the same size.
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSize';
        msg = sprintf('Inputs must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end