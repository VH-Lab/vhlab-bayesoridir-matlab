% DiffFromUniform.m
% This function should be placed in a directory structure like:
% .../+vis/+bayes/+assess/DiffFromUniform.m

function totalAbsoluteDifference = DiffFromUniform(values, probabilities)
%DiffFromUniform Computes the integrated absolute difference between a PDF and a uniform PDF.
%
%   totalAbsDiff = vis.bayes.assess.DiffFromUniform(values, probabilities)
%   Calculates the total integrated absolute difference between a piece-wise
%   constant Probability Density Function (PDF) derived from the input
%   'values' and 'probabilities', and a uniform PDF defined over the same
%   'values'.
%
%   The difference is quantified as:
%       integral | PDF_actual(x) - PDF_uniform(x) | dx
%
%   The function uses 'vis.bayes.assess.pieceWiseConstantPDF' to determine
%   the bin edges, widths, and PDF heights for both the actual distribution
%   and the reference uniform distribution. The total absolute difference
%   is then computed by summing the contributions from each bin:
%       Sum_i [ |height_actual_i - height_uniform_i| * width_i ]
%
%   Inputs:
%       values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
%                     mustBeStrictlyIncreasingForDiff(values)}
%           A column vector of unique values representing the points around which
%           PDF bins are centered. Must be strictly increasing. Passed to
%           vis.bayes.assess.pieceWiseConstantPDF.
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSizeForDiff(values, probabilities)}
%           A column vector of non-negative probabilities (or unnormalized
%           masses/weights) corresponding to each element in 'values'. Passed to
%           vis.bayes.assess.pieceWiseConstantPDF.
%
%   Outputs:
%       totalAbsoluteDifference (1,1) double
%           The calculated total integrated absolute difference.
%           - A value of 0 indicates the actual PDF is identical to the uniform PDF.
%           - Larger positive values indicate greater deviation from uniformity.
%           - The maximum possible value is 2 (when the distributions have no overlap
%             and both integrate to 1).
%           - Returns NaN if inputs are invalid (e.g., empty values/probabilities
%             leading to an undefined PDF).
%
%   Requires:
%       - vis.bayes.assess.pieceWiseConstantPDF.m
%
%   Example 1: Actual distribution is uniform
%       vals = [10; 20; 30];
%       probsActual = [1; 1; 1]; % Will be normalized to [1/3; 1/3; 1/3]
%       diff = vis.bayes.assess.DiffFromUniform(vals, probsActual);
%       % Expected diff = 0 (or very close due to floating point)
%
%   Example 2: Actual distribution is a spike, very different from uniform
%       vals = [10; 20; 30];
%       probsActual = [0.01; 0.98; 0.01]; % Highly concentrated at 20
%       diff = vis.bayes.assess.DiffFromUniform(vals, probsActual);
%       % Expected diff to be a positive value, approaching 2.
%       % PDF_actual heights: [0.001, 0.098, 0.001] (for widths 10, if edges are 5,15,25,35)
%       % PDF_uniform heights: [1/30, 1/30, 1/30] approx [0.0333, 0.0333, 0.0333]
%       % Bin1: |0.001 - 0.0333|*10 = |-0.0323|*10 = 0.323
%       % Bin2: |0.098 - 0.0333|*10 = |0.0647|*10  = 0.647
%       % Bin3: |0.001 - 0.0333|*10 = |-0.0323|*10 = 0.323
%       % Total diff approx 0.323 + 0.647 + 0.323 = 1.293
%
%   Example 3: Single point distribution
%       vals = [10];
%       probsActual = [1];
%       diff = vis.bayes.assess.DiffFromUniform(vals, probsActual);
%       % For a single point, the actual and uniform PDFs are identical.
%       % Expected diff = 0

    arguments
        values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
                     mustBeStrictlyIncreasingForDiff(values)}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSizeForDiff(values, probabilities)}
    end

    totalAbsoluteDifference = NaN; % Default output for invalid scenarios

    % --- Step 1: Get PDF model for the actual probabilities ---
    try
        [~, binWidths, ~, pdfHeightsActual] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);
    catch ME
        error('vis:bayes:assess:DiffFromUniform:ActualPDFError', ...
              'Error modeling PDF for actual probabilities: %s.', ME.message);
    end

    % Validate outputs from pieceWiseConstantPDF for actual distribution
    if isempty(pdfHeightsActual) || isempty(binWidths) || (numel(pdfHeightsActual) ~= numel(binWidths))
        error('vis:bayes:assess:DiffFromUniform:InvalidActualPDFModel', ...
              'PDF model from pieceWiseConstantPDF for actual probabilities is invalid or empty.');
    end

    % --- Step 2: Create uniform probabilities and get PDF model for uniform distribution ---
    numValues = numel(values);
    if numValues == 0
        % This should be caught by mustBeStrictlyIncreasingForDiff if it enforces non-empty.
        % If 'values' is empty, pieceWiseConstantPDF would have warned/errored or returned empty.
        error('vis:bayes:assess:DiffFromUniform:NoValues', ...
              'Input "values" is empty. Cannot compute difference from uniform.');
    end
    
    uniformProbabilities = ones(numValues, 1); % Unnormalized, pieceWiseConstantPDF will normalize

    try
        % For the uniform distribution, binWidths and binEdges will be the same
        % as for the actual distribution since 'values' are the same.
        [~, ~, ~, pdfHeightsUniform] = vis.bayes.assess.pieceWiseConstantPDF(values, uniformProbabilities);
    catch ME
        error('vis:bayes:assess:DiffFromUniform:UniformPDFError', ...
              'Error modeling PDF for uniform probabilities: %s.', ME.message);
    end
    
    % Validate outputs from pieceWiseConstantPDF for uniform distribution
    if isempty(pdfHeightsUniform) || (numel(pdfHeightsUniform) ~= numValues)
         error('vis:bayes:assess:DiffFromUniform:InvalidUniformPDFModel', ...
              'PDF model from pieceWiseConstantPDF for uniform probabilities is invalid or empty.');
    end

    % --- Step 3: Calculate Total Integrated Absolute Difference ---
    % Sum_i [ |height_actual_i - height_uniform_i| * width_i ]
    
    % Check if either PDF resulted in all zero heights (e.g. sum of input probabilities was zero)
    % If actual probabilities summed to zero, pdfHeightsActual will be all zeros.
    % If uniform probabilities summed to zero (not possible with ones()), pdfHeightsUniform would be zeros.
    % The comparison still makes sense. If one is all zeros and the other is not,
    % the difference will be sum(abs(nonZeroHeights) .* widths).
    % If both are all zeros, difference is 0.

    differencesInHeight = abs(pdfHeightsActual - pdfHeightsUniform);
    totalAbsoluteDifference = sum(differencesInHeight .* binWidths);

end

% Local validation functions
function mustBeStrictlyIncreasingForDiff(v)
    if isempty(v) && numel(v)==0
        % Allow empty to pass this validator; main function will error if numValues is 0.
        return;
    end
    if ~(isvector(v) && issorted(v, 'strictascend'))
        if numel(v) == 1 && isnumeric(v) % Single element is okay
            return;
        end
        eid = 'validation:mustBeStrictlyIncreasingForDiff';
        msg = 'Input "values" for DiffFromUniform calculation must be a strictly increasing vector.';
        throwAsCaller(MException(eid, msg));
    end
end

function mustHaveSameSizeForDiff(a, b)
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSizeForDiff';
        msg = sprintf('Inputs "values" and "probabilities" for DiffFromUniform calculation must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end
