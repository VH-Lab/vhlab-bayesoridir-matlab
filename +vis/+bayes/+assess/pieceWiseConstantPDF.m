% pieceWiseConstantPDF.m
% This function should be placed in a directory structure like:
% .../+vis/+bayes/+assess/pieceWiseConstantPDF.m

function [normalizedBinMasses, binWidths, binEdges, pdfHeights] = pieceWiseConstantPDF(values, probabilities)
%pieceWiseConstantPDF Constructs a piece-wise constant PDF model.
%
%   [normMasses, widths, edges, heights] = pieceWiseConstantPDF(values, probabilities)
%   Models a probability distribution as a piece-wise constant Probability
%   Density Function (PDF) based on a series of discrete values and their
%   associated probabilities (or unnormalized weights).
%
%   The function assumes each input value in 'values' represents the conceptual
%   center of a bin. The probability mass associated with this value is
%   distributed uniformly across this bin.
%
%   Bin Definition:
%   The bin associated with each 'values(i)' is defined by extending halfway
%   to the preceding value 'values(i-1)' and halfway to the succeeding value
%   'values(i+1)'.
%   - For an interior value 'values(i)', its bin spans from
%     (values(i-1) + values(i))/2 to (values(i) + values(i+1))/2.
%   - For the first value 'values(1)', its bin starts at
%     values(1) - (values(2) - values(1))/2 and ends at
%     (values(1) + values(2))/2.
%   - For the last value 'values(end)', its bin starts at
%     (values(end-1) + values(end))/2 and ends at
%     values(end) + (values(end) - values(end-1))/2.
%   - If only a single value is provided (i.e., numel(values) == 1), its bin
%     is centered at this value with a default width of 1.0.
%
%   Normalization:
%   The input 'probabilities' are treated as masses or weights for each bin.
%   These are normalized so that the sum of 'normalizedBinMasses' (the total
%   area under the PDF) equals 1.0, unless the sum of input probabilities
%   is zero or negligible.
%
%   Inputs:
%       values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
%                     mustBeStrictlyIncreasing(values)}
%           A column vector of unique values representing the points around which
%           PDF bins are centered. Must be strictly increasing.
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSize(values, probabilities)}
%           A column vector of non-negative probabilities (or unnormalized
%           masses/weights) corresponding to each element in 'values'.
%
%   Outputs:
%       normalizedBinMasses (:,1) double
%           A column vector representing the normalized probability mass within
%           each bin. sum(normalizedBinMasses) will be 1.0 if the input
%           probabilities sum to a positive value. Corresponds element-wise
%           to the input 'values'.
%
%       binWidths (:,1) double
%           A column vector of the calculated widths for each bin. Corresponds
%           element-wise to the input 'values'.
%
%       binEdges (N+1,1) double
%           A column vector of bin edges, where N is the number of input values
%           (numel(values)). binEdges(i) is the left edge of the bin associated
%           with values(i), and binEdges(i+1) is its right edge. Thus, there
%           are N+1 edges for N bins.
%
%       pdfHeights (:,1) double
%           A column vector representing the height of the PDF within each bin.
%           Calculated as normalizedBinMasses ./ binWidths. Corresponds
%           element-wise to the input 'values'.
%
%   Example 1: Standard case with three points
%       vals = [10; 20; 30];
%       probs = [1; 2; 1]; % Unnormalized weights, sum = 4
%       [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(vals, probs);
%       % Expected normMasses: [0.25; 0.5; 0.25]
%       % Expected edges: [5; 15; 25; 35]
%       % Expected widths: [10; 10; 10]
%       % Expected heights: [0.025; 0.05; 0.025]
%
%   Example 2: Single point input
%       valsSingle = [10];
%       probsSingle = [5];
%       [normMassesS, widthsS, edgesS, heightsS] = vis.bayes.assess.pieceWiseConstantPDF(valsSingle, probsSingle);
%       % Expected normMassesS: [1]
%       % Expected edgesS: [9.5; 10.5]
%       % Expected widthsS: [1]
%       % Expected heightsS: [1]

    arguments
        values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
                     mustBeStrictlyIncreasing(values)}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSize(values, probabilities)}
    end

    % Initialize outputs for early return scenarios
    normalizedBinMasses = double.empty(0,1);
    binWidths = double.empty(0,1);
    binEdges = double.empty(0,1);
    pdfHeights = double.empty(0,1);
    tolerance = sqrt(eps); % Tolerance for floating point comparisons

    % --- Input Handling and Preprocessing ---
    if isempty(values)
        warning('MATLAB:pieceWiseConstantPDF:EmptyInput', 'Input values vector is empty.');
        return;
    end

    % Ensure inputs are column vectors (already enforced by (:,1) but good for clarity)
    values = values(:);
    probabilities = probabilities(:);
    numValues = numel(values);

    % --- Normalize Probabilities (to get Bin Masses) ---
    totalProbabilityMass = sum(probabilities);
    if totalProbabilityMass <= tolerance
        warning('pieceWiseConstantPDF:ZeroTotalProbability', ...
                'Total sum of input probabilities is zero or negligible. Normalized masses and PDF heights will be zero.');
        normalizedBinMasses = zeros(numValues, 1);
    else
        normalizedBinMasses = probabilities / totalProbabilityMass;
    end

    % --- Calculate Bin Edges and Widths ---
    if numValues == 1
        % Special handling for a single value: bin centered with width 1 by convention
        binEdges = [values(1) - 0.5; values(1) + 0.5];
        binWidths = [1.0];
    else % numValues >= 2
        midpoints = (values(1:end-1) + values(2:end)) / 2.0;
        firstBinLeftEdge = values(1) - (values(2) - values(1)) / 2.0;
        lastBinRightEdge = values(end) + (values(end) - values(end-1)) / 2.0;
        binEdges = [firstBinLeftEdge; midpoints; lastBinRightEdge];
        binWidths = diff(binEdges);

        if any(binWidths <= tolerance)
            warning('pieceWiseConstantPDF:ZeroWidthBin', ...
                    'Bin calculation resulted in zero or negligible widths. This may occur if input values are distinct but extremely close. PDF heights might be affected.');
        end
    end

    % --- Calculate PDF Heights ---
    pdfHeights = zeros(numValues, 1); % Initialize
    if totalProbabilityMass > tolerance % Only calculate heights if there's mass
        % Calculate heights for bins with non-negligible width
        validWidthIndices = binWidths > tolerance;
        if any(validWidthIndices)
            pdfHeights(validWidthIndices) = normalizedBinMasses(validWidthIndices) ./ binWidths(validWidthIndices);
        end
    end
    
    % Ensure all outputs are column vectors
    normalizedBinMasses = normalizedBinMasses(:);
    binWidths = binWidths(:);
    binEdges = binEdges(:);
    pdfHeights = pdfHeights(:);

end

% Local validation function for 'values' input
function mustBeStrictlyIncreasing(v)
    if isempty(v) % Allow empty to pass, main function handles it
        return;
    end
    if ~(isvector(v) && issorted(v, 'strictascend'))
        if numel(v) == 1 && isnumeric(v) % Single element is trivially strictly increasing
            return;
        end
        eid = 'validation:mustBeStrictlyIncreasing';
        msg = 'Input "values" must be a strictly increasing vector.';
        throwAsCaller(MException(eid, msg));
    end
end

% Local validation function (required by arguments block for probabilities)
function mustHaveSameSize(a, b)
%MUSTHAVESAMESIZE Validate that inputs have the same size.
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSize';
        msg = sprintf('Inputs must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end
