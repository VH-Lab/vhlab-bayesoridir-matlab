% expectedSquaredError.m
% .../+vis/+bayes/+assess/expectedSquaredError.m

function ese = expectedSquaredError(values, probabilities, trueValue)
%expectedSquaredError Computes the expected squared error for a piece-wise constant PDF.
%
%   ese = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue)
%   Calculates the expected squared error (ESE) between a random variable X,
%   whose distribution is modeled by a piece-wise constant PDF derived from
%   'values' and 'probabilities', and a given 'trueValue'.
%
%   The ESE is defined as E[(X - trueValue)^2] = integral( (x - trueValue)^2 * PDF(x) dx ).
%
%   This function uses 'vis.bayes.assess.pieceWiseConstantPDF' to determine
%   the bin edges, widths, and PDF heights of the modeled distribution.
%   The integral is then computed by summing the contributions from each bin.
%
%   Inputs:
%       values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
%                     mustBeStrictlyIncreasingForESE(values)}
%           A column vector of unique values representing the points around which
%           PDF bins are centered. Must be strictly increasing. Passed to
%           vis.bayes.assess.pieceWiseConstantPDF.
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSizeForESE(values, probabilities)}
%           A column vector of non-negative probabilities (or unnormalized
%           masses/weights) corresponding to each element in 'values'. Passed to
%           vis.bayes.assess.pieceWiseConstantPDF.
%
%       trueValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
%           A scalar representing the true value against which the squared
%           error is calculated.
%
%   Outputs:
%       ese (1,1) double
%           The calculated expected squared error. Returns NaN if inputs are
%           invalid (e.g., empty values/probabilities leading to an undefined PDF)
%           or if the sum of probabilities is zero.
%
%   Requires:
%       - vis.bayes.assess.pieceWiseConstantPDF.m
%       (And its associated local validation functions if not on path, though this
%        function defines its own for clarity or if they differ slightly).
%
%   Example 1: Simple symmetric distribution, trueValue at center
%       vals = [10; 20; 30];
%       probs = [1; 2; 1]; % PDF heights [0.025; 0.05; 0.025] on bins of width 10
%       tv = 20;
%       err = vis.bayes.assess.expectedSquaredError(vals, probs, tv);
%       % err should be the variance of this distribution.
%       % For bin 1 (5 to 15, h=0.025): (0.025/3)*((15-20)^3 - (5-20)^3) = (0.025/3)*(-125 - (-3375)) = (0.025/3)*3250 = 27.0833
%       % For bin 2 (15 to 25, h=0.05): (0.05/3)*((25-20)^3 - (15-20)^3) = (0.05/3)*(125 - (-125)) = (0.05/3)*250 = 4.1667
%       % For bin 3 (25 to 35, h=0.025): (0.025/3)*((35-20)^3 - (25-20)^3) = (0.025/3)*(3375 - 125) = (0.025/3)*3250 = 27.0833
%       % Total ESE = 27.0833 + 4.1667 + 27.0833 = 58.3333
%
%   Example 2: Single point distribution
%       vals = [10];
%       probs = [1];
%       tv = 10;
%       err = vis.bayes.assess.expectedSquaredError(vals, probs, tv); % ESE should be small, related to bin width
%       % Bin is [9.5, 10.5], height 1.
%       % (1/3)*((10.5-10)^3 - (9.5-10)^3) = (1/3)*(0.5^3 - (-0.5)^3) = (1/3)*(0.125 - (-0.125)) = (1/3)*0.25 = 0.08333
%
%       tv_offset = 11;
%       err_offset = vis.bayes.assess.expectedSquaredError(vals, probs, tv_offset);
%       % (1/3)*((10.5-11)^3 - (9.5-11)^3) = (1/3)*((-0.5)^3 - (-1.5)^3) = (1/3)*(-0.125 - (-3.375)) = (1/3)*3.25 = 1.08333

    arguments
        values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
                     mustBeStrictlyIncreasingForESE(values)}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSizeForESE(values, probabilities)}
        trueValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    ese = NaN; % Default output for invalid scenarios

    % --- Step 1: Model the PDF using pieceWiseConstantPDF ---
    try
        % Assuming vis.bayes.assess.pieceWiseConstantPDF is on the path
        [~, ~, binEdges, pdfHeights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);
    catch ME
        warning('vis:bayes:assess:expectedSquaredError:PDFHelperError', ...
                'Error calling pieceWiseConstantPDF: %s. Cannot compute ESE.', ME.message);
        return; % ESE remains NaN
    end

    % --- Step 2: Validate results from pieceWiseConstantPDF ---
    if isempty(pdfHeights) || isempty(binEdges) || numel(binEdges) ~= (numel(pdfHeights) + 1)
        % This can happen if pieceWiseConstantPDF returned empty due to empty inputs,
        % or if there's an unexpected structure.
        % pieceWiseConstantPDF should warn for empty inputs itself.
        warning('vis:bayes:assess:expectedSquaredError:InvalidPDFModel', ...
                'PDF model from pieceWiseConstantPDF is invalid or empty. Cannot compute ESE.');
        return; % ESE remains NaN
    end
    
    % If all pdfHeights are zero (e.g., sum of input probabilities was zero), ESE is 0.
    if all(abs(pdfHeights) < sqrt(eps))
        ese = 0;
        return;
    end

    % --- Step 3: Calculate Expected Squared Error ---
    numBins = numel(pdfHeights);
    binEseContributions = zeros(numBins, 1);

    for i = 1:numBins
        li = binEdges(i);     % Left edge of the current bin
        ri = binEdges(i+1);   % Right edge of the current bin
        hi = pdfHeights(i);   % PDF height in the current bin

        if hi < sqrt(eps) % Skip bins with (near) zero probability height
            binEseContributions(i) = 0;
            continue;
        end
        
        % Integral of (x - trueValue)^2 * hi dx over [li, ri]
        % = hi * [ (x - trueValue)^3 / 3 ] from li to ri
        term_ri = (ri - trueValue)^3;
        term_li = (li - trueValue)^3;
        
        integral_val = (hi / 3.0) * (term_ri - term_li);
        binEseContributions(i) = integral_val;
    end

    ese = sum(binEseContributions);

end

% Local validation functions (can be identical to those in pieceWiseConstantPDF)
% It's good practice for a function to define its own validators or ensure
% they are accessible from a shared utility location.
function mustBeStrictlyIncreasingForESE(v)
    if isempty(v) % Allow empty to pass to pieceWiseConstantPDF, which handles it
        return;
    end
    if ~(isvector(v) && issorted(v, 'strictascend'))
        if numel(v) == 1 && isnumeric(v) % Single element is trivially strictly increasing
            return;
        end
        eid = 'validation:mustBeStrictlyIncreasingForESE';
        msg = 'Input "values" for ESE calculation must be a strictly increasing vector.';
        throwAsCaller(MException(eid, msg));
    end
end

function mustHaveSameSizeForESE(a, b)
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSizeForESE';
        msg = sprintf('Inputs "values" and "probabilities" for ESE calculation must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end
