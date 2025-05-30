% TestExpectedSquaredError.m
% This test file should be located such that the package +vis can be accessed.
% For example, if expectedSquaredError.m is in Project/+vis/+bayes/+assess/expectedSquaredError.m,
% this test file could be in Project/+vis/+bayes/+assess/TestExpectedSquaredError.m

classdef TestExpectedSquaredError < matlab.unittest.TestCase
    % Test class for vis.bayes.assess.expectedSquaredError

    properties
        AbsTol = 1e-9; % Absolute tolerance for floating point comparisons
        RelTol = 1e-6; % Relative tolerance for some comparisons
    end

    methods (Test)

        function testSymmetricDistributionCenterTrueValue(testCase)
            % Test case 1: From function example 1
            % Symmetric distribution, trueValue at the center of symmetry.
            % ESE should be the variance of the distribution.
            values = [10; 20; 30];
            probabilities = [1; 2; 1]; % Sum = 4
            trueValue = 20;
            
            % Manual calculation from doc:
            % Bin 1 (5 to 15, h=0.025): (0.025/3)*((15-20)^3 - (5-20)^3) = (0.025/3)*3250 = 27.08333333
            % Bin 2 (15 to 25, h=0.05): (0.05/3)*((25-20)^3 - (15-20)^3) = (0.05/3)*250  = 4.16666667
            % Bin 3 (25 to 35, h=0.025): (0.025/3)*((35-20)^3 - (25-20)^3) = (0.025/3)*3250 = 27.08333333
            % Total ESE = 27.08333333 + 4.16666667 + 27.08333333 = 58.33333333
            expectedESE = 58.33333333333333;

            actualESE = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
            testCase.verifyEqual(actualESE, expectedESE, 'AbsTol', testCase.AbsTol, ...
                'ESE for symmetric distribution with trueValue at center failed.');
        end

        function testSymmetricDistributionOffsetTrueValue(testCase)
            % Test case 2: Symmetric distribution, trueValue offset from center.
            values = [10; 20; 30];
            probabilities = [1; 2; 1]; % Sum = 4
            trueValue = 15; % Offset true value
            
            % PDF heights: [0.025, 0.05, 0.025] for bins [5,15], [15,25], [25,35]
            % Bin 1 (5 to 15, h=0.025), T=15: (0.025/3)*((15-15)^3 - (5-15)^3) = (0.025/3)*(0 - (-1000)) = (0.025/3)*1000 = 8.33333333
            % Bin 2 (15 to 25, h=0.05), T=15: (0.05/3)*((25-15)^3 - (15-15)^3) = (0.05/3)*(1000 - 0) = (0.05/3)*1000 = 16.66666667
            % Bin 3 (25 to 35, h=0.025), T=15: (0.025/3)*((35-15)^3 - (25-15)^3) = (0.025/3)*(20^3 - 10^3) = (0.025/3)*(8000 - 1000) = (0.025/3)*7000 = 58.33333333
            % Total ESE = 8.33333333 + 16.66666667 + 58.33333333 = 83.33333333
            expectedESE = 83.33333333333333;

            actualESE = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
            testCase.verifyEqual(actualESE, expectedESE, 'AbsTol', testCase.AbsTol, ...
                'ESE for symmetric distribution with offset trueValue failed.');
        end

        function testSinglePointAtTrueValue(testCase)
            % Test case 3: From function example 2, part 1
            % Single point distribution, trueValue at the point.
            values = [10];
            probabilities = [1];
            trueValue = 10;
            
            % Bin is [9.5, 10.5], height 1.
            % (1/3)*((10.5-10)^3 - (9.5-10)^3) = (1/3)*(0.5^3 - (-0.5)^3) = (1/3)*(0.125 + 0.125) = (1/3)*0.25
            expectedESE = 0.25 / 3.0; % 0.083333333333333

            actualESE = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
            testCase.verifyEqual(actualESE, expectedESE, 'AbsTol', testCase.AbsTol, ...
                'ESE for single point at trueValue failed.');
        end

        function testSinglePointOffsetTrueValue(testCase)
            % Test case 4: From function example 2, part 2
            % Single point distribution, trueValue offset from the point.
            values = [10];
            probabilities = [1];
            trueValue = 11;
            
            % Bin is [9.5, 10.5], height 1.
            % (1/3)*((10.5-11)^3 - (9.5-11)^3) = (1/3)*((-0.5)^3 - (-1.5)^3) = (1/3)*(-0.125 - (-3.375)) = (1/3)*3.25
            expectedESE = 3.25 / 3.0; % 1.083333333333333

            actualESE = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
            testCase.verifyEqual(actualESE, expectedESE, 'AbsTol', testCase.AbsTol, ...
                'ESE for single point with offset trueValue failed.');
        end

        function testZeroProbabilitiesInput(testCase)
            % Test case 5: Probabilities are all zero.
            values = [10; 20; 30];
            probabilities = [0; 0; 0];
            trueValue = 15;
            expectedESE = 0; % If no probability mass, ESE should be 0.

            % pieceWiseConstantPDF will warn about zero total probability.
            % expectedSquaredError should handle this by returning ESE = 0.
            testCase.verifyWarning(@() vis.bayes.assess.pieceWiseConstantPDF(values, probabilities), ...
                                   'pieceWiseConstantPDF:ZeroTotalProbability');
            
            actualESE = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
            testCase.verifyEqual(actualESE, expectedESE, 'AbsTol', testCase.AbsTol, ...
                'ESE for zero probabilities input failed.');
        end

        function testEmptyInputs(testCase)
            % Test case 6: Empty values and probabilities.
            valuesEmpty = double.empty(0,1);
            probabilitiesEmpty = double.empty(0,1);
            trueValue = 10;
            % expectedESE = NaN; % Function returns NaN for this case

            % pieceWiseConstantPDF warns for empty input.
            % expectedSquaredError then warns about invalid PDF model and returns NaN.
            testCase.verifyWarning(@() vis.bayes.assess.expectedSquaredError(valuesEmpty, probabilitiesEmpty, trueValue), ...
                                   'vis:bayes:assess:expectedSquaredError:InvalidPDFModel');
            
            actualESE = vis.bayes.assess.expectedSquaredError(valuesEmpty, probabilitiesEmpty, trueValue);
            testCase.verifyTrue(isnan(actualESE), 'ESE should be NaN for empty inputs.');
        end
        
        function testUniformDistribution(testCase)
            % Test case 7: A simple uniform distribution.
            % PDF = 1/10 for x in [0, 10], 0 otherwise.
            % Let's model this with two points.
            values = [2.5; 7.5]; % Centers of two halves
            probabilities = [0.5; 0.5]; % Equal mass
            trueValue = 5; % Mean of the uniform distribution
            
            % pieceWiseConstantPDF for values=[2.5; 7.5], probs=[0.5;0.5]
            % Edges: L(2.5)-(7.5-2.5)/2 = 2.5-2.5=0
            %        M(2.5+7.5)/2 = 5
            %        R(7.5)+(7.5-2.5)/2 = 7.5+2.5=10
            % So, edges = [0; 5; 10].
            % Widths = [5; 5].
            % NormMasses = [0.5; 0.5].
            % PDFHeights = [0.5/5; 0.5/5] = [0.1; 0.1]. This is PDF = 1/10 on [0,10].

            % ESE = Var(X) for trueValue = E[X]
            % For U(a,b), Var(X) = (b-a)^2 / 12. Here U(0,10), Var(X) = (10-0)^2/12 = 100/12 = 25/3.
            expectedESE = 25.0/3.0; % 8.333333333333333

            actualESE = vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue);
            testCase.verifyEqual(actualESE, expectedESE, 'AbsTol', testCase.AbsTol, ...
                'ESE for uniform distribution failed.');
        end

        % --- Tests for Input Validation and Error Handling ---
        function testValuesNotStrictlyIncreasingError(testCase)
            values = [10; 5]; % Not strictly increasing
            probabilities = [1; 1];
            trueValue = 7;
            testCase.verifyError(@() vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue), ...
                'validation:mustBeStrictlyIncreasingForESE');
        end

        function testProbabilitiesSizeMismatchError(testCase)
            values = [10; 20];
            probabilities = [1]; % Size mismatch
            trueValue = 15;
            testCase.verifyError(@() vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue), ...
                'validation:mustHaveSameSizeForESE');
        end
        
        function testNegativeProbabilitiesErrorPropagates(testCase)
            % This error is caught by pieceWiseConstantPDF, which expectedSquaredError calls.
            % The error from pieceWiseConstantPDF should propagate.
            values = [10; 20];
            probabilities = [1; -1]; % Negative probability
            trueValue = 15;
            testCase.verifyError(@() vis.bayes.assess.expectedSquaredError(values, probabilities, trueValue), ...
                'MATLAB:validators:mustBeNonnegative');
        end
        
        function testNonFiniteTrueValueError(testCase)
            values = [10; 20];
            probabilities = [1; 1];
            
            trueValueNaN = NaN;
            testCase.verifyError(@() vis.bayes.assess.expectedSquaredError(values, probabilities, trueValueNaN), ...
                'MATLAB:validators:mustBeFinite', 'Test with NaN trueValue.');
            
            trueValueInf = Inf;
            testCase.verifyError(@() vis.bayes.assess.expectedSquaredError(values, probabilities, trueValueInf), ...
                'MATLAB:validators:mustBeFinite', 'Test with Inf trueValue.');
        end
        
        function testValuesNotVectorErrorPropagates(testCase)
            % Error from pieceWiseConstantPDF due to shape validation in its arguments block
            valuesMatrix = [10 20; 30 40];
            probabilitiesMatrix = [1 1; 1 1]; % Matching size for its own validation
            trueValue = 25;
            
            % The error 'MATLAB:validation:IncompatibleSize' comes from pieceWiseConstantPDF's
            % 'values (:,1)' argument validation.
            testCase.verifyError(@() vis.bayes.assess.expectedSquaredError(valuesMatrix, probabilitiesMatrix, trueValue), ...
                'MATLAB:validation:IncompatibleSize');
        end

    end
end
