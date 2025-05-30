function inPercentile = isInCentralPercentile2(values, probabilities, testValue, testPercentiles, options)
%isInCentralPercentile Determines if a value is within the central N% of a distribution.
%
%   inPercentile = isInCentralPercentile2(values, probabilities, testValue, testPercentiles)
%   Determines if testValue lies within the central N percent of a probability
%   distribution defined by 'values' and 'probabilities'. The "center" of the
%   distribution is defined by its mode. The distribution is modeled as a
%   piece-wise constant PDF based on the input points.
%
%   inPercentile = isInCentralPercentile2(..., Name, Value) specifies options
%   using one or more Name, Value pair arguments.
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
%           The scalar value to test against the central percentile intervals.
%
%       testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                             mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
%           A vector of percentages (0 to 100).
%
%   Optional Name-Value Pairs:
%       PlotType ('none' | 'pdf' | 'cdf')
%           'none' (default): No plot.
%           'pdf': Plots the PDF, mode, interval boundaries, and testValue.
%           'cdf': Plots the CDF, mode, interval boundaries, and testValue.
%
%   Outputs:
%       inPercentile (:,1) logical
%           True if testValue falls within the central N% interval.
%
%   Central Interval Definition (Revised):
%   The interval starts with the bin containing the mode. It is then expanded
%   iteratively by adding adjacent bins (alternating sides when possible).
%   If expansion on one side hits the edge of the distribution's support,
%   expansion continues on the other side until the total probability mass
%   within the interval reaches N/100. If adding a full bin exceeds N/100,
%   only the necessary fraction of that bin (assuming uniform density within it)
%   is added to reach exactly N/100.
%
%   Example:
%       values = [10; 20; 30]; % Mode is 10
%       probabilities = [0.6; 0.3; 0.1];
%       testVal = 18;
%       percentiles = [50; 70; 95];
%       % For 70%, starts with bin 1 (0.6). Needs 0.1 more. Adds fraction of bin 2.
%       isInside = isInCentralPercentile(values, probabilities, testVal, percentiles, 'PlotType', 'pdf');
%       % Expect isInside = [false; true; true]

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
    tol = sqrt(eps); % Tolerance for floating point comparisons


    % --- Input Handling and Preprocessing ---
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
        modeValue = values(1);
        totalProb = sum(probabilities);
         if totalProb <= tol % Check if probability is essentially zero
             warning('Total probability mass is zero or negative.');
             inPercentile = false(nPerc, 1);
             return;
         end
        for i = 1:nPerc
           if testPercentiles(i) > tol
               lower_bounds(i) = modeValue;
               upper_bounds(i) = modeValue;
               inPercentile(i) = (abs(testValue - modeValue) < tol); % Check with tolerance
           else
               lower_bounds(i) = NaN;
               upper_bounds(i) = NaN;
               inPercentile(i) = false;
           end
        end
        if ~strcmpi(options.PlotType, 'none')
            warning('Plotting is not well-defined for a single-point distribution. Skipping plot.');
        end
        return;
    end

    % Normalize probabilities to sum to 1
    totalProb = sum(probabilities);
    if totalProb <= tol
        warning('Total probability mass is zero or negative.');
        inPercentile = false(nPerc, 1);
        return;
    elseif abs(totalProb - 1.0) > tol
        normalized_probabilities = probabilities / totalProb;
    else
        normalized_probabilities = probabilities;
    end

    % --- Define PDF Bins ---
    midpoints = (values(1:end-1) + values(2:end)) / 2;
    first_half_width = (midpoints(1) - values(1));
    last_half_width = (values(end) - midpoints(end));
    first_edge = values(1);
    last_edge = values(end);
    if nVal > 1 && first_edge >= midpoints(1) - tol
        first_edge = values(1) - eps(values(1));
    end
     if nVal > 1 && last_edge <= midpoints(end) + tol
         last_edge = values(end) + eps(values(end));
     end
    edges = [first_edge; midpoints; last_edge];
    widths = diff(edges);
    if any(widths <= tol)
        bad_width_indices = find(widths <= tol);
        problematic_value_indices = unique([bad_width_indices; bad_width_indices+1]);
        problematic_value_indices(problematic_value_indices > nVal) = [];
        warning('Bin calculation resulted in zero or negative widths near values: %s. Check input values.', num2str(values(problematic_value_indices)'));
        inPercentile = false(nPerc, 1);
        return;
    end

    % --- Calculate PDF and CDF (CDF needed for plotting option) ---
    pdf_heights = normalized_probabilities ./ widths;
    cdf_vals = zeros(nVal + 1, 1); % Still needed for CDF plot option
    cdf_vals(2:end) = cumsum(normalized_probabilities);
    cdf_vals(end) = 1.0;

    % --- Find Mode ---
    [~, mode_idx] = max(normalized_probabilities);
    mode_idx = find(abs(normalized_probabilities - normalized_probabilities(mode_idx)) < tol, 1, 'first'); % Handle ties
    modeValue = values(mode_idx);

    % --- Calculate Intervals using Growing Window Method ---
    for i = 1:nPerc
        percentile = testPercentiles(i);
        target_mass = percentile / 100.0;

        if target_mass <= tol
            lower_bounds(i) = NaN;
            upper_bounds(i) = NaN;
            inPercentile(i) = false;
            continue;
        elseif target_mass >= 1.0 - tol
             lower_bounds(i) = edges(1);
             upper_bounds(i) = edges(end);
             % Check inclusion below
        else
            % --- Start Growing Window Method ---
            prob_mode_bin = normalized_probabilities(mode_idx);
            pdf_mode_bin = pdf_heights(mode_idx);

            if prob_mode_bin >= target_mass - tol
                % Mode bin contains enough mass. Find sub-interval within mode bin.
                if pdf_mode_bin > tol
                    sub_width = target_mass / pdf_mode_bin;
                    % Center around modeValue (center of the bin value)
                    lb = modeValue - sub_width / 2;
                    ub = modeValue + sub_width / 2;
                    % Clamp to mode bin edges
                    lower_bounds(i) = max(edges(mode_idx), lb);
                    upper_bounds(i) = min(edges(mode_idx + 1), ub);
                else % Mode bin has mass but zero PDF height (e.g., infinite width?) -> problematic case
                     warning('isInCentralPercentile:ZeroPDFHeight',...
                             'Mode bin %d has mass %.4g but zero PDF height. Assigning full bin edges.', mode_idx, prob_mode_bin);
                     lower_bounds(i) = edges(mode_idx);
                     upper_bounds(i) = edges(mode_idx+1);
                end
            else
                % Start expansion from mode bin
                current_mass = prob_mode_bin;
                lower_bin_idx = mode_idx;
                upper_bin_idx = mode_idx;
                current_lower_bound = edges(mode_idx);
                current_upper_bound = edges(mode_idx + 1);

                % Alternating expansion (start lower)
                side_to_expand = -1; % -1 for lower, +1 for upper

                while current_mass < target_mass - tol % Loop until target mass reached

                    try_lower = (side_to_expand == -1);
                    try_upper = (side_to_expand == 1);

                    can_go_lower = (lower_bin_idx > 1);
                    can_go_upper = (upper_bin_idx < nVal);

                    made_progress = false;

                    % Determine which sides are *possible* to expand to
                    expand_lower_possible = can_go_lower;
                    expand_upper_possible = can_go_upper;

                    % Determine which side to *attempt* based on alternation & possibility
                    attempt_lower = try_lower && expand_lower_possible;
                    attempt_upper = try_upper && expand_upper_possible;

                    % If preferred side is not possible, try the other side
                    if try_lower && ~attempt_lower && expand_upper_possible
                        attempt_upper = true;
                        attempt_lower = false; % Don't try lower again this iteration
                    elseif try_upper && ~attempt_upper && expand_lower_possible
                        attempt_lower = true;
                        attempt_upper = false; % Don't try upper again this iteration
                    end

                    % --- Try Expanding Lower ---
                    if attempt_lower
                        next_lower = lower_bin_idx - 1;
                        prob_to_add = normalized_probabilities(next_lower);
                        if current_mass + prob_to_add >= target_mass - tol
                            needed_mass = target_mass - current_mass;
                            if prob_to_add > tol
                                fraction = needed_mass / prob_to_add;
                                current_lower_bound = edges(next_lower + 1) - fraction * widths(next_lower);
                            else % Adding zero-prob bin needed_mass must be zero
                                 current_lower_bound = edges(next_lower); % take full zero bin? unlikely needed
                            end
                            current_mass = target_mass; % Reached target
                            made_progress = true;
                        else
                            current_mass = current_mass + prob_to_add;
                            current_lower_bound = edges(next_lower);
                            lower_bin_idx = next_lower;
                            made_progress = true;
                        end
                        side_to_expand = 1; % Switch preferred side
                        if current_mass >= target_mass - tol, break; end % Check if done
                    end

                    % --- Try Expanding Upper ---
                    % Re-evaluate attempt possibility if lower attempt failed but upper is possible
                     if ~made_progress && try_lower && ~attempt_lower && expand_upper_possible
                          attempt_upper = true; % Force attempt upper if lower failed/impossible
                     end

                    if attempt_upper
                        next_upper = upper_bin_idx + 1;
                        prob_to_add = normalized_probabilities(next_upper);
                         if current_mass + prob_to_add >= target_mass - tol
                            needed_mass = target_mass - current_mass;
                            if prob_to_add > tol
                                 fraction = needed_mass / prob_to_add;
                                 current_upper_bound = edges(next_upper) + fraction * widths(next_upper);
                            else
                                 current_upper_bound = edges(next_upper+1);
                            end
                            current_mass = target_mass; % Reached target
                            made_progress = true;
                        else
                            current_mass = current_mass + prob_to_add;
                            current_upper_bound = edges(next_upper + 1);
                            upper_bin_idx = next_upper;
                            made_progress = true;
                        end
                        side_to_expand = -1; % Switch preferred side
                        if current_mass >= target_mass - tol, break; end % Check if done
                    end

                    % --- Check for stall ---
                    if ~made_progress && ~expand_lower_possible && ~expand_upper_possible
                        % Cannot expand further, but haven't reached target mass
                        warning('isInCentralPercentile:ExpansionStall', ...
                                'Could not reach target mass %.4f. Expansion stalled at %.4f.', target_mass, current_mass);
                        break; % Exit loop, use current bounds
                    end

                end % end while loop

                lower_bounds(i) = current_lower_bound;
                upper_bounds(i) = current_upper_bound;
            end % end if/else for mode bin check / expansion
        end % end if/else for target_mass range

        % Check if testValue is within the calculated interval
        if ~isnan(lower_bounds(i)) && ~isnan(upper_bounds(i))
            % Add tolerance to avoid issues exactly at the boundary
             inPercentile(i) = (testValue >= lower_bounds(i) - tol) && (testValue <= upper_bounds(i) + tol);
        else
            inPercentile(i) = false;
        end
    end % End of loop i = 1:nPerc


    % --- Optional Plotting (No changes needed here) ---
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
             if maxY < tol, maxY = 1; end % Handle zero heights case
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
        valid_lb = lower_bounds(~isnan(lower_bounds));
        valid_ub = upper_bounds(~isnan(upper_bounds));
        unique_bounds = unique([valid_lb; valid_ub]);

        if ~isempty(unique_bounds)
             bound_style = '--';
             bound_color = [0.6 0.6 0.6];
             xline(unique_bounds, bound_style, 'Color', bound_color, 'LineWidth', 1, 'HandleVisibility', 'off');
             h_bound_leg = xline(NaN, bound_style, 'Color', bound_color, 'LineWidth', 1); % Invisible line for legend
             plotHandles(end+1) = h_bound_leg;
             plotLabels{end+1} = 'Central Percentile Bounds';
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
        legend(plotHandles, plotLabels, 'Location', 'best');
        x_range = edges(end) - edges(1);
         if abs(x_range) < tol
             x_range = abs(edges(1)) * 0.2 + 1;
         end
        xlim_low = edges(1) - 0.05 * x_range;
        xlim_high = edges(end) + 0.05 * x_range;
        all_vals_of_interest = [testValue; modeValue; unique_bounds]; % Include bounds in auto-ranging
        xlim_low = min(xlim_low, min(all_vals_of_interest) - 0.05 * x_range);
        xlim_high = max(xlim_high, max(all_vals_of_interest) + 0.05 * x_range);
        xlim([xlim_low xlim_high]);

        if strcmpi(options.PlotType, 'pdf')
             ylim([0 maxY]);
        else % CDF
             ylim([-0.05 1.05]);
        end
        hold off;
    end % End plotting block

end % End main function

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

