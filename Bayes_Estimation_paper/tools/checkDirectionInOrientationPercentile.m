function inPercentile = checkDirectionInOrientationPercentile(direction_angles, probabilities, testDirection, testPercentiles, options)
%checkDirectionInOrientationPercentile Checks if a direction is in the central N% of its derived orientation distribution.
%
%   inPercentile = checkDirectionInOrientationPercentile(direction_angles, probabilities, testDirection, testPercentiles, ...)
%   Takes a distribution defined over direction angles [0, 360), transforms
%   it into a distribution over orientation angles [0, 180) which is then
%   rotated so the peak probability aligns with angle 0. It then checks if
%   the provided testDirection (also transformed into the rotated orientation
%   space) falls within the central N% of this final distribution for each
%   N specified in testPercentiles.
%
%   This function utilizes two helper functions:
%   1. transformDirectionToOrientation: To convert the direction data to the
%      rotated orientation space and find the transformed test direction.
%   2. isInCentralPercentile2: To determine if the transformed test direction
%      lies within the central N% of the rotated orientation distribution.
%
%   Inputs:
%       direction_angles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
%                     mustBeInRange(direction_angles, 0, 360, 'exclude-upper')}
%           A vector of original direction angles in degrees [0, 360).
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSize(direction_angles, probabilities)}
%           A vector of probabilities corresponding to each direction angle.
%
%       testDirection (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                            mustBeInRange(testDirection, 0, 360, 'exclude-upper')}
%           The specific direction angle (degrees, [0, 360)) to test.
%
%       testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                             mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
%           A vector of percentages (0 to 100) defining the central intervals
%           to check against in the rotated orientation distribution.
%
%   Optional Name-Value Pairs:
%       PlotDirectional (1,1) logical = false
%           If true, generates the polar plot and rotated orientation plot
%           from the 'transformDirectionToOrientation' function call.
%
%       PlotOrientationPercentile (1,1) logical = false
%           If true, generates the PDF plot showing the central percentile
%           bounds from the 'isInCentralPercentile2' function call relative
%           to the rotated orientation distribution.
%
%   Outputs:
%       inPercentile (:,1) logical
%           A logical vector, the same size as 'testPercentiles'. Each element
%           is true if the transformed testDirection falls within the central
%           interval defined by the corresponding percentage in testPercentiles
%           for the rotated orientation distribution. Returns false if the
%           transformation or percentile check could not be performed.
%
%   Requires:
%       - transformDirectionToOrientation.m
%       - isInCentralPercentile2.m
%       (And their associated helper/validation functions if not defined locally)
%
%   Examples:
%       % Example 1: Simple Bimodal Directional Data
%       angles1 = [45; 135; 225; 315]; % Directions with symmetry
%       probs1  = [0.1; 0.4; 0.1; 0.4];  % Probabilities -> Peak Orientation 135 deg
%       testDir1 = 45;                  % Maps to orientation 45 deg
%       percentiles1 = [10; 50; 90];    % Check multiple percentile bounds
%
%       % The peak orientation is 135 deg. Rotated distribution has values at
%       % 0 deg (prob 0.8 from 135/315) and -90 deg (prob 0.2 from 45/225).
%       % Test direction 45 maps to orientation 45, rotated to -90 deg.
%       % Check if -90 is within the central N% of this {-90:0.2, 0:0.8} distribution.
%       isInside1 = checkDirectionInOrientationPercentile(angles1, probs1, testDir1, percentiles1, ...
%              'PlotDirectional', true, 'PlotOrientationPercentile', true);
%       disp('Example 1 Results:');
%       disp(table(percentiles1, isInside1));
%       % Expected: isInside1 = [false; false; true] (approx.)
%
%       % Example 2: Different Peak Location
%       angles2 = [10; 190; 90; 270];   % Directions
%       probs2  = [0.4; 0.4; 0.1; 0.1];  % Probabilities -> Peak Orientation 10 deg
%       testDir2 = 90;                  % Maps to orientation 90 deg
%       percentiles2 = [50; 95];
%
%       % Peak orientation is 10 deg. Rotated distribution is at 0 deg (prob 0.8 from 10/190)
%       % and 80 deg (prob 0.2 from 90/270).
%       % Test direction 90 maps to orientation 90, rotated angle 80.
%       % Check if 80 is within the central N% of this {0:0.8, 80:0.2} distribution.
%       isInside2 = checkDirectionInOrientationPercentile(angles2, probs2, testDir2, percentiles2, ...
%              'PlotDirectional', true, 'PlotOrientationPercentile', true);
%       disp('Example 2 Results:');
%       disp(table(percentiles2, isInside2));
%       % Expected: isInside2 = [false; true] (approx.)
%
%       % Example 3: More complex distribution with plotting
%       angles3 = (0:30:330)'; % Directions every 30 deg
%       % Assign probabilities (higher near 0/180 and 90/270)
%       p3 = [0.15; 0.1; 0.05; 0.05; 0.05; 0.1; 0.15; 0.1; 0.05; 0.05; 0.05; 0.1];
%       probs3 = p3 / sum(p3); % Normalize
%       testDir3 = 100;           % Test direction near 90/270 peaks
%       percentiles3 = [80];
%
%       % Check percentile and generate both plots
%       isInside3 = checkDirectionInOrientationPercentile(angles3, probs3, testDir3, percentiles3, ...
%                     'PlotDirectional', true, 'PlotOrientationPercentile', true);
%       disp('Example 3 Results (Percentile=80):');
%       disp(isInside3); % Display true/false result
%       % Two figures should be generated for visual inspection.

    arguments
        direction_angles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
                              mustBeInRange(direction_angles, 0, 360, 'exclude-upper')}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                               mustBeVector, mustBeNonnegative, mustHaveSameSize(direction_angles, probabilities)}
        testDirection (1,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                               mustBeInRange(testDirection, 0, 360, 'exclude-upper')}
        testPercentiles (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                                 mustBeVector, mustBeInRange(testPercentiles, 0, 100)}
        options.PlotDirectional (1,1) logical = false
        options.PlotOrientationPercentile (1,1) logical = false
    end

    % --- Step 1: Transform Direction Data to Rotated Orientation Space ---

    % Call the first function, passing the testDirection and plotting option
    try
        % Ensure transformDirectionToOrientation is on the path or defined below
        [mod_angles, mod_probs, ~, trans_dir_val] = transformDirectionToOrientation(...
            direction_angles, probabilities, ...
            'directionValue', testDirection, ...
            'MakePlot', options.PlotDirectional);
    catch ME
        warning('Error during call to transformDirectionToOrientation: %s', ME.message);
        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
             error('Function transformDirectionToOrientation.m not found. Ensure it is on the MATLAB path.');
        end
        inPercentile = false(size(testPercentiles)); % Return default failure
        return;
    end

    % --- Step 2: Validate Transformation Results ---

    % Check if transformation yielded valid results needed for percentile check
    if isempty(mod_angles) || isempty(mod_probs) || isnan(trans_dir_val)
        warning('Transformation resulted in empty data or invalid transformed test direction. Cannot check percentile.');
        inPercentile = false(size(testPercentiles)); % Return default failure
        return;
    end

    % Check if probabilities sum to a non-zero value (isInCentralPercentile2 also checks this)
    if sum(mod_probs) <= sqrt(eps)
         warning('Summed orientation probabilities are zero or negligible. Cannot check percentile.');
         inPercentile = false(size(testPercentiles));
         return;
    end

    % Handle the case where only one unique orientation resulted
    if numel(mod_angles) == 1
         fprintf('Note: Only one unique orientation angle resulted from transformation. Percentile check might be trivial.\n');
    end


    % --- Step 3: Check Percentile Inclusion in Rotated Orientation Space ---

    % Determine plot type for the second function call
    plotType = 'none';
    if options.PlotOrientationPercentile
        plotType = 'pdf'; % Defaulting to PDF plot for percentile check
    end

    % Call the second function using the transformed data
    try
        % Ensure isInCentralPercentile2 is on the path or defined below
        inPercentile = isInCentralPercentile2(...
            mod_angles, mod_probs, ...    % Rotated orientation angles and their probabilities
            trans_dir_val, ...            % The transformed test value
            testPercentiles, ...          % The percentiles to check
            'PlotType', plotType);        % Pass plot option
    catch ME
         warning('Error during call to isInCentralPercentile2: %s', ME.message);
          if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
             error('Function isInCentralPercentile2.m not found. Ensure it is on the MATLAB path.');
          end
         inPercentile = false(size(testPercentiles)); % Return default failure
         return;
    end

end % End main function

% --------------------------------------------------------------------------
% Local Helper / Validation Functions (if not on path)
% --------------------------------------------------------------------------
% NOTE: The following validation function is needed by the arguments block.
% Ensure it is defined here OR accessible on the MATLAB path (e.g., if it
% was defined locally within transformDirectionToOrientation.m, it needs to be
% defined here again or made globally accessible).
% Similarly for mustHaveSameSize.

function mustHaveSameSize(a, b)
%MUSTHAVESAMESIZE Validate that inputs have the same size.
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSize';
        msg = sprintf('Inputs must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end

% function mustBeValidDirectionOrNaN(val)
% %MUSTBEVALIDDIRECTIONORNAN Validate value is NaN or in [0, 360).
%     if ~(isnan(val) || (val >= 0 && val < 360))
%         eid = 'validation:mustBeValidDirectionOrNaN';
%         msg = 'Input must be NaN or a direction angle in the range [0, 360).';
%         throwAsCaller(MException(eid, msg));
%     end
% end

% Assumes transformDirectionToOrientation and isInCentralPercentile2 are M-files
% accessible on the MATLAB path.