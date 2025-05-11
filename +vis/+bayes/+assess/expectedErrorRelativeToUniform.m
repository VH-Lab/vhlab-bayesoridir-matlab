% expectedErrorRelativeToUniform.m
% .../+vis/+bayes/+assess/expectedErrorRelativeToUniform.m

function relErr = expectedErrorRelativeToUniform(values, probabilities, trueValue)
%expectedErrorRelativeToUniform Computes ESE relative to a uniform distribution.
%
%   relErr = vis.bayes.assess.expectedErrorRelativeToUniform(values, probabilities, trueValue)
%   Calculates the expected squared error (ESE) of a given probability
%   distribution and compares it to the ESE of a uniform distribution
%   defined over the same 'values'. The result is the relative difference:
%   (ESE_actual - ESE_uniform) / ESE_uniform.
%
%   The function uses 'vis.bayes.assess.expectedSquaredError' as a helper
%   to compute both the ESE for the provided 'probabilities' (ESE_actual)
%   and the ESE for a uniform probability distribution over the 'values'
%   (ESE_uniform).
%
%   Inputs:
%       values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
%                     mustBeStrictlyIncreasingForRelESE(values)}
%           A column vector of unique values representing the points around which
%           PDF bins are centered. Must be strictly increasing.
%
%       probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
%                           mustBeVector, mustBeNonnegative, mustHaveSameSizeForRelESE(values, probabilities)}
%           A column vector of non-negative probabilities (or unnormalized
%           masses/weights) corresponding to each element in 'values'.
%
%       trueValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
%           A scalar representing the true value against which the squared
%           error is calculated for both distributions.
%
%   Outputs:
%       relErr (1,1) double
%           The relative error: (ESE_actual - ESE_uniform) / ESE_uniform.
%           - Returns 0 if ESE_uniform is near zero and ESE_actual is also near zero.
%           - Throws an error if ESE_uniform is near zero but ESE_actual is not.
%           - Throws an error if either ESE calculation fails or results in NaN.
%
%   Requires:
%       - vis.bayes.assess.expectedSquaredError.m
%       - vis.bayes.assess.pieceWiseConstantPDF.m (as a dependency of expectedSquaredError)
%
%   Example 1: Actual distribution is better than uniform
%       vals = [10; 20; 30];
%       probsActual = [0.1; 0.8; 0.1]; tv = 20;
%       relErr1 = vis.bayes.assess.expectedErrorRelativeToUniform(vals, probsActual, tv);
%       % Expected relErr1 approx -0.8846
%
%   Example 2: Actual distribution is worse than uniform
%       vals = [10; 20; 30];
%       probsActual = [0.4; 0.2; 0.4]; tv = 20;
%       relErr2 = vis.bayes.assess.expectedErrorRelativeToUniform(vals, probsActual, tv);
%       % Expected relErr2 approx 0.4615
%
%   Example 3: Single point (Actual and Uniform are effectively the same)
%       vals = [10]; probsActual = [1]; tv = 10;
%       relErr3 = vis.bayes.assess.expectedErrorRelativeToUniform(vals, probsActual, tv);
%       % Expected relErr3 = 0

    arguments
        values (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector, ...
                     mustBeStrictlyIncreasingForRelESE(values)}
        probabilities (:,1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                           mustBeVector, mustBeNonnegative, mustHaveSameSizeForRelESE(values, probabilities)}
        trueValue (1,1) {mustBeNumeric, mustBeReal, mustBeFinite}
    end

    % relErr = NaN; % No longer default to NaN, will error out instead.
    tolerance = sqrt(eps); % For checking if ESE_uniform is near zero

    % --- Step 1: Calculate ESE for the actual given probabilities ---
    try
        eseActual = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
    catch ME
        error('vis:bayes:assess:expectedErrorRelativeToUniform:ActualESECalcError', ...
              'Error calculating ESE for actual probabilities: %s. Cannot compute relative error.', ME.message);
    end

    if isnan(eseActual)
        error('vis:bayes:assess:expectedErrorRelativeToUniform:ActualESENaN', ...
              'ESE for actual probabilities resulted in NaN. Cannot compute relative error.');
    end

    % --- Step 2: Create uniform probabilities and calculate ESE_uniform ---
    numValues = numel(values);
    if numValues == 0
        % This should ideally be caught by the validator for 'values' if it enforces non-empty.
        % If 'values' can be empty and pass validation, this error is necessary.
        error('vis:bayes:assess:expectedErrorRelativeToUniform:NoValuesForUniform', ...
              'Input "values" is empty. Cannot create uniform distribution or compute relative error.');
    end

    uniformProbabilities = ones(numValues, 1) / numValues; % Normalized uniform probabilities

    try
        eseUniform = vis.bayes.assess.expectedSquaredError(values, uniformProbabilities, trueValue);
    catch ME
        error('vis:bayes:assess:expectedErrorRelativeToUniform:UniformESECalcError', ...
              'Error calculating ESE for uniform probabilities: %s. Cannot compute relative error.', ME.message);
    end
    
    if isnan(eseUniform)
        error('vis:bayes:assess:expectedErrorRelativeToUniform:UniformESENaN', ...
              'ESE for uniform probabilities resulted in NaN. Cannot compute relative error.');
    end

    % --- Step 3: Calculate Relative Error ---
    if abs(eseUniform) < tolerance
        % Denominator (ESE_uniform) is close to zero
        if abs(eseActual) < tolerance
            % Numerator (ESE_actual) also close to zero (0/0 case) -> define as 0
            relErr = 0;
        else
            % Numerator non-zero, denominator zero -> undefined/infinite relative error
            error('vis:bayes:assess:expectedErrorRelativeToUniform:ZeroUniformESE', ...
                  'ESE for uniform distribution is near zero (%.4g) while actual ESE (%.4g) is not, leading to undefined or infinite relative error.', ...
                  eseUniform, eseActual);
        end
    else
        % Standard calculation
        relErr = (eseActual - eseUniform) / eseUniform;
    end

end

% Local validation functions
function mustBeStrictlyIncreasingForRelESE(v)
    % Validator for 'values'. It should ensure 'v' is suitable for numel(v) > 0 if used later.
    % If 'values' must be non-empty, this validator should enforce it.
    % For this function, numel(values) == 0 is explicitly checked and errored.
    if isempty(v) && numel(v)==0 
        % Allow empty to pass this validator; main function will error if numValues is 0.
        % This is to avoid erroring here if 'values' is allowed to be empty by some calling contexts
        % before this specific function's logic (numValues == 0 check) is hit.
        % However, if 'values' must always be non-empty for this function,
        % this validator should throw an error for isempty(v).
        % Given the explicit numValues == 0 check, this is okay.
        return;
    end
    if ~(isvector(v) && issorted(v, 'strictascend'))
        if numel(v) == 1 && isnumeric(v) % Single element is okay
            return;
        end
        eid = 'validation:mustBeStrictlyIncreasingForRelESE';
        msg = 'Input "values" for relative ESE calculation must be a strictly increasing vector.';
        throwAsCaller(MException(eid, msg));
    end
end

function mustHaveSameSizeForRelESE(a, b)
    if ~isequal(size(a), size(b))
        eid = 'validation:mustHaveSameSizeForRelESE';
        msg = sprintf('Inputs "values" and "probabilities" for relative ESE calculation must have the same size. Size A = [%s], Size B = [%s].', ...
                      num2str(size(a)), num2str(size(b)));
        throwAsCaller(MException(eid, msg));
    end
end
