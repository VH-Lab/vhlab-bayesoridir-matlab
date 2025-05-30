% TestDiffFromUniform.m
% This test file should be located such that the package +vis can be accessed.
% For example, if DiffFromUniform.m is in Project/+vis/+bayes/+assess/DiffFromUniform.m,
% this test file could be in Project/+vis/+bayes/+assess/TestDiffFromUniform.m

classdef TestDiffFromUniform < matlab.unittest.TestCase
    % Test class for vis.bayes.assess.DiffFromUniform

    properties
        AbsTol = 1e-9; % Absolute tolerance for floating point comparisons
        % RelTol = 1e-6; % Relative tolerance, if needed
    end

    methods (Test)

        function testActualIsUniform(testCase)
            % Test case 1: Actual distribution is already uniform.
            % From function example 1.
            values = [10; 20; 30];
            probabilities = [1; 1; 1]; % These will be normalized to uniform by pieceWiseConstantPDF
            
            expectedDifference = 0.0;

            actualDifference = vis.bayes.assess.DiffFromUniform(values, probabilities);
            testCase.verifyEqual(actualDifference, expectedDifference, 'AbsTol', testCase.AbsTol, ...
                'Difference should be zero when actual distribution is uniform.');
        end

        function testHighlyConcentratedActual(testCase)
            % Test case 2: Actual distribution is a spike, very different from uniform.
            % From function example 2.
            values = [10; 20; 30];
            probabilitiesActual = [0.01; 0.98; 0.01]; % Sum = 1, highly concentrated at 20
            
            expectedDifference = 1.2933333333333333;

            actualDifference = vis.bayes.assess.DiffFromUniform(values, probabilitiesActual);
            testCase.verifyEqual(actualDifference, expectedDifference, 'AbsTol', testCase.AbsTol, ...
                'Difference for highly concentrated actual distribution is incorrect.');
        end

        function testSinglePointDistribution(testCase)
            % Test case 3: Single point distribution.
            % From function example 3.
            % For a single point, the actual and uniform PDFs derived by pieceWiseConstantPDF are identical.
            values = [10];
            probabilitiesActual = [1]; % Any non-zero probability
            
            expectedDifference = 0.0;

            actualDifference = vis.bayes.assess.DiffFromUniform(values, probabilitiesActual);
            testCase.verifyEqual(actualDifference, expectedDifference, 'AbsTol', testCase.AbsTol, ...
                'Difference for a single point distribution should be zero.');
        end

        function testZeroActualProbabilities(testCase)
            % Test case 4: Actual probabilities are all zero.
            % pdfHeightsActual will be all zeros.
            % pdfHeightsUniform will be standard uniform heights (non-zero if numValues > 0).
            % Difference should be sum(|0 - height_uniform_i| * width_i) = sum(height_uniform_i * width_i) = 1.
            values = [10; 20; 30];
            probabilitiesActual = [0; 0; 0];
            
            expectedDifference = 1.0;

            % pieceWiseConstantPDF for actual will warn about zero total probability.
            testCase.verifyWarning(@() vis.bayes.assess.pieceWiseConstantPDF(values, probabilitiesActual), ...
                                   'pieceWiseConstantPDF:ZeroTotalProbability');

            actualDifference = vis.bayes.assess.DiffFromUniform(values, probabilitiesActual);
            testCase.verifyEqual(actualDifference, expectedDifference, 'AbsTol', testCase.AbsTol, ...
                'Difference when actual probabilities are zero is incorrect.');
        end

        function testBimodalActualDistribution(testCase)
            % Test case 5: A bimodal distribution compared to uniform.
            values = [10; 20; 30; 40; 50];
            probabilitiesActual = [0.4; 0.05; 0.1; 0.05; 0.4]; % Peaks at 10 and 50
            
            actualDifference = vis.bayes.assess.DiffFromUniform(values, probabilitiesActual);
            testCase.verifyGreaterThan(actualDifference, 0, ...
                'Difference for bimodal distribution should be positive.');
            testCase.verifyLessThanOrEqual(actualDifference, 2.0, ...
                'Difference should not exceed the theoretical maximum of 2.');
            
            expectedDifferenceBimodal = 0.8;
            testCase.verifyEqual(actualDifference, expectedDifferenceBimodal, 'AbsTol', testCase.AbsTol);
        end

        % --- Tests for Input Validation and Error Handling ---
        function testEmptyValuesError(testCase)
            % Test case 6: 'values' is empty.
            valuesEmpty = double.empty(0,1);
            probabilities = double.empty(0,1); % Must match size for other validator
            
            % When values is empty, pieceWiseConstantPDF returns empty pdfHeightsActual.
            % This leads to the 'InvalidActualPDFModel' error in DiffFromUniform.
            testCase.verifyError(@() vis.bayes.assess.DiffFromUniform(valuesEmpty, probabilities), ...
                'vis:bayes:assess:DiffFromUniform:InvalidActualPDFModel', ...
                'Error for empty values input not thrown or incorrect ID.');
        end
        
        function testValuesNotStrictlyIncreasingError(testCase)
            % Test case 7: 'values' not strictly increasing.
            values = [10; 5]; 
            probabilities = [1; 1]; % Must match size
            testCase.verifyError(@() vis.bayes.assess.DiffFromUniform(values, probabilities), ...
                'validation:mustBeStrictlyIncreasingForDiff');
        end

        function testProbabilitiesSizeMismatchError(testCase)
            % Test case 8: 'values' and 'probabilities' size mismatch.
            values = [10; 20];
            probabilities = [1]; 
            testCase.verifyError(@() vis.bayes.assess.DiffFromUniform(values, probabilities), ...
                'validation:mustHaveSameSizeForDiff');
        end
        
        function testNegativeProbabilitiesErrorPropagates(testCase)
            % Test case 9: Negative probabilities.
            % This is caught by the 'mustBeNonnegative' validator in DiffFromUniform's arguments block.
            values = [10; 20];
            probabilities = [1; -1]; 
            
            testCase.verifyError(@() vis.bayes.assess.DiffFromUniform(values, probabilities), ...
                'MATLAB:validators:mustBeNonnegative');
        end
        
        function testNonFiniteValuesErrorPropagates(testCase)
            % Test case 10: Non-finite values (e.g., NaN, Inf)
            % This is caught by the 'mustBeFinite' validator for 'values' in DiffFromUniform's arguments block.
            valuesNaN = [10; NaN; 30];
            probabilities = [1; 1; 1]; % Must match size
            testCase.verifyError(@() vis.bayes.assess.DiffFromUniform(valuesNaN, probabilities), ...
                'MATLAB:validators:mustBeFinite');

            valuesInf = [10; Inf; 30];
             probabilitiesInf = [1;1;1]; % Must match size
            testCase.verifyError(@() vis.bayes.assess.DiffFromUniform(valuesInf, probabilitiesInf), ...
                'MATLAB:validators:mustBeFinite');
        end
    end
end
