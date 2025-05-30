function [modified_angles, modified_probabilities, orientation_angles, transformedDirectionValue] = transformDirectionToOrientation(angles, probabilities, options)
%transformDirectionToOrientation Convert direction dist to rotated orientation dist.
%
%   [mod_angles, mod_probs, orient_angles, trans_dir] = transformDirectionToOrientation(angles, probabilities, ...)
%   Transforms a probability distribution defined over direction angles
%   (0 <= angles < 360 degrees) to a distribution over orientation angles
%   (0 <= angles < 180 degrees). Also optionally transforms a specific
%   direction angle to the new rotated orientation space.
%
%   The transformation involves:
%   1. Mapping each direction angle to its equivalent orientation angle using
%      mod(angle, 180).
%   2. Summing the probabilities for all direction angles that map to the
%      same orientation angle.
%   3. Identifying the orientation angle(s) with the maximum summed probability.
%      If multiple peaks exist, the one with the smallest angle is chosen.
%   4. Rotating the resulting orientation angles so the chosen peak orientation
%      corresponds to zero. The rotated angles represent the difference from
%      the peak orientation, wrapped into the range [-90, 90) degrees.
%   5. If 'directionValue' is provided, it is mapped to its equivalent orientation
%      and then rotated relative to the peak orientation using the same method.
%
%   Inputs:
%       angles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
%                     mustBeInRange(angles, 0, 360, 'exclude-upper')}
%           A vector of direction angles in degrees. Values must be in the
%           range [0, 360).
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSize(angles, probabilities)}
%           A vector of probabilities corresponding to each element in 'angles'.
%           Must have the same size as 'angles' and contain non-negative values.
%
%   Optional Name-Value Pairs:
%       MakePlot (1,1) logical = false
%           Specifies whether to generate a figure with two subplots:
%           - Panel 1: Polar plot of the original direction angles and probabilities.
%           - Panel 2: Stem plot of the rotated orientation angles vs. summed probabilities.
%           If 'directionValue' is also provided, it will be indicated on the plots.
%
%       directionValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                             mustBeValidDirectionOrNaN(options.directionValue)} = NaN
%           A specific direction angle (in the range [0, 360)) to be transformed
%           into the rotated orientation space. If provided (i.e., not NaN),
%           the fourth output argument 'transformedDirectionValue' will be calculated.
%           Defaults to NaN, indicating no specific value is provided.
%
%   Outputs:
%       modified_angles (:,1) double
%           The rotated orientation angles in degrees, representing the difference
%           from the peak orientation, wrapped into the range [-90, 90). Size
%           corresponds to the number of unique orientation angles.
%
%       modified_probabilities (:,1) double
%           The summed probabilities corresponding to each element in
%           modified_angles (and orientation_angles).
%
%       orientation_angles (:,1) double
%           The unique orientation angles (degrees, [0, 180)) before rotation,
%           corresponding element-wise to modified_probabilities.
%
%       transformedDirectionValue (1,1) double
%           The rotated orientation angle corresponding to the input 'directionValue'.
%           Calculated only if 'directionValue' is provided (not NaN).
%           Returns NaN if 'directionValue' was not provided or if the transformation
%           could not be completed (e.g., empty input).
%
%   Example 1: Simple peak, track a direction
%       angles = [45; 135; 225; 315]; % Directions
%       probs = [0.1; 0.4; 0.1; 0.4];  % Probabilities
%       [mod_a, mod_p, orient_a, trans_dir] = transformDirectionToOrientation(angles, probs, 'directionValue', 225);
%       % orient_a = [45; 135], mod_p = [0.2; 0.8], Peak=135, mod_a = [-90; 0]
%       % Direction 225 maps to orientation 45.
%       % trans_dir will be -90 (45 deg orientation is -90 deg from peak 135).
%
%   Example 2: Uniform orientation, plotting enabled, track a direction
%       angles = (0:45:315)';       % Directions every 45 deg
%       probs = repmat(0.125, 8, 1); % Uniform probability
%       [~, ~, ~, trans_dir] = transformDirectionToOrientation(angles, probs, 'MakePlot', true, 'directionValue', 180);
%       % Peak orientation is 0 deg.
%       % Direction 180 maps to orientation 0.
%       % trans_dir will be 0.
%       % A figure with two panels will be displayed, showing the tracked direction.

    arguments
        angles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
                     mustBeInRange(angles, 0, 360, 'exclude-upper')}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSize(angles, probabilities)}
        options.MakePlot (1,1) logical = false
        options.directionValue (1,1) {mustBeNumeric, mustBeReal, ...
                                      mustBeValidDirectionOrNaN(options.directionValue)} = NaN
    end

    transformedDirectionValue = NaN; % Initialize fourth output

    % --- Handle Empty Input ---
    if isempty(angles)
        modified_angles = [];
        modified_probabilities = [];
        orientation_angles = [];
        % transformedDirectionValue remains NaN
        if options.MakePlot
            warning('Input angles are empty. Cannot generate plot.');
        end
        return;
    end

    tol = sqrt(eps); % Tolerance for floating point comparisons
    directionValueProvided = ~isnan(options.directionValue); % Check if optional arg was given

    % --- Step 1: Map Directions to Orientations ---
    orient_equiv = mod(angles, 180);

    % --- Step 2: Sum Probabilities for Each Orientation ---
    [unique_orientations_sorted, ~, ic] = unique(orient_equiv);
    summed_probabilities = accumarray(ic, probabilities, [], @sum);
    orientation_angles = unique_orientations_sorted; % Assign output

    if ~iscolumn(orientation_angles), orientation_angles = orientation_angles'; end
    if ~iscolumn(summed_probabilities), summed_probabilities = summed_probabilities'; end

    if isempty(orientation_angles) || isempty(summed_probabilities)
        modified_angles = [];
        modified_probabilities = [];
        % transformedDirectionValue remains NaN
        if options.MakePlot
            warning('Could not determine valid orientation angles or probabilities. Cannot generate plot.');
        end
        return;
    end

    % --- Step 3: Find the Peak Orientation ---
    max_prob = max(summed_probabilities);
    peak_indices = find(abs(summed_probabilities - max_prob) < tol * max(1, max_prob));

    if isempty(peak_indices)
         warning('Could not find peak probability. Using first orientation angle as reference.');
         peak_orientation = orientation_angles(1);
    else
        % Choose the first peak index (corresponding to the smallest angle among peaks)
        peak_orientation = orientation_angles(peak_indices(1));
    end

    % --- Step 4: Rotate Orientations relative to Peak ---
    angle_diff = orientation_angles - peak_orientation;
    modified_angles = mod(angle_diff + 90, 180) - 90;
    modified_probabilities = summed_probabilities; % Assign output

    % --- Step 5: Transform optional directionValue ---
    if directionValueProvided
        directionValue = options.directionValue;
        orient_equiv_dir = mod(directionValue, 180);
        angle_diff_dir = orient_equiv_dir - peak_orientation;
        transformedDirectionValue = mod(angle_diff_dir + 90, 180) - 90; % Assign output
    end

    % --- Step 6: Optional Plotting ---
    if options.MakePlot
        figure;
        plotHandles = []; % Store handles for legend
        plotLabels = {};  % Store labels for legend

        % Panel 1: Original Directions (Polar Plot)
        subplot(1, 2, 1);
        % Use user's modification to close the loop
        h_data = polarplot(deg2rad([angles(:); angles(1)]), [probabilities(:); probabilities(1)], ...
                           'o-', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
        plotHandles(end+1) = h_data(1); % Use first handle from plot object array if necessary
        plotLabels{end+1} = 'Data';
        hold on; % Hold on to add direction value ray

        % Add tracked direction if provided
        if directionValueProvided
             max_r_disp = max(probabilities);
             if max_r_disp < tol, max_r_disp = 1; end
             h_dir = polarplot([0 deg2rad(directionValue)], [0 max_r_disp], 'g--', 'LineWidth', 1.5);
             plotHandles(end+1) = h_dir;
             plotLabels{end+1} = sprintf('Direction (%.1f^\\circ)', directionValue);
        end

        ax = gca;
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        max_r_lim = max(probabilities);
        if max_r_lim < tol, max_r_lim = 1; end
        rlim([0, max_r_lim * 1.1]);
        title('Original Directional Distribution');
        legend(plotHandles, plotLabels, 'Location', 'bestoutside');
        hold off;

        % Clear handles for next plot
        plotHandles = [];
        plotLabels = {};

        % Panel 2: Rotated Orientations (Stem Plot)
        subplot(1, 2, 2);
        hold on; % Hold on to add lines
        h_stem = stem(modified_angles, modified_probabilities, 'filled', 'LineWidth', 1.5, 'Color', 'b');
        plotHandles(end+1) = h_stem;
        plotLabels{end+1} = 'Orient. Probs';

        % Add vertical line for transformed direction if provided
        if directionValueProvided
            h_trans_dir = xline(transformedDirectionValue, 'g--', 'LineWidth', 1.5);
            plotHandles(end+1) = h_trans_dir;
            plotLabels{end+1} = sprintf('Transformed Dir. (%.1f^\\circ)', transformedDirectionValue);
        end

        grid on;
        xlabel('Angle Difference from Peak Orientation (degrees)');
        ylabel('Summed Probability');
        title(sprintf('Rotated Orientation Dist. (Peak at %.1f^\\circ)', peak_orientation));
        max_abs_angle = max(abs(modified_angles));
        if isempty(max_abs_angle), max_abs_angle = 0; end % Handle empty case
        lim_val = max(max_abs_angle, 90);
        if directionValueProvided % Ensure tracked value is visible
            lim_val = max(lim_val, abs(transformedDirectionValue));
        end
        lim = lim_val * 1.1;
        xlim([-lim, lim]);
        max_y = max(modified_probabilities);
        if isempty(max_y) || max_y < tol, max_y = 1; end
        ylim([0 max_y*1.1]);
        legend(plotHandles, plotLabels, 'Location', 'best');
        hold off;
    end

end % End main function

% --------------------------------------------------------------------------
% Local Helper Functions
% --------------------------------------------------------------------------

% Custom validation function for directionValue (allows NaN or range)
function mustBeValidDirectionOrNaN(val)
%MUSTBEVALIDDIRECTIONORNAN Validate value is NaN or in [0, 360).
    if ~(isnan(val) || (val >= 0 && val < 360))
        eid = 'validation:mustBeValidDirectionOrNaN';
        msg = 'Input must be NaN or a direction angle in the range [0, 360).';
        throwAsCaller(MException(eid, msg));
    end
end

% Custom validation function (required by arguments block)
function mustHaveSameSize(a, b)
%MUSTHAVESAMESIZE Validate that inputs have the same size.
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSize';
        msg = sprintf('Inputs must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end