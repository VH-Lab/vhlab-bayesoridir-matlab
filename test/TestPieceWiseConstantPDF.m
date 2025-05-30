% TestPieceWiseConstantPDF.m
% This test file should be located such that the package +vis can be accessed.
% For example, if pieceWiseConstantPDF is in Project/+vis/+bayes/+assess/pieceWiseConstantPDF.m,
% this test file could be in Project/+vis/+bayes/+assess/TestPieceWiseConstantPDF.m

classdef TestPieceWiseConstantPDF < matlab.unittest.TestCase
    % Test class for vis.bayes.assess.pieceWiseConstantPDF

    properties
        AbsTol = 1e-9; % Absolute tolerance for floating point comparisons
        RelTol = 1e-6; % Relative tolerance for some comparisons (e.g., for very large numbers)
    end

    methods (Test)

        function testSingleValue(testCase)
            % Test case 1: Single value input
            values = [20.0];
            probabilities = [10.0]; % Unnormalized
            
            expectedNormMasses = [1.0];
            expectedBinEdges = [19.5; 20.5]; % Default width 1, centered
            expectedBinWidths = [1.0];
            expectedPdfHeights = [1.0];

            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges, expectedBinEdges, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths, expectedBinWidths, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(heights, expectedPdfHeights, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(sum(heights .* widths), 1.0, 'AbsTol', testCase.AbsTol, ...
                'Total PDF area should integrate to 1.0.');
        end

        function testTwoValues(testCase)
            % Test case 2: Two values input
            values = [5.0; 15.0];
            probabilities = [0.3; 0.7]; % Already normalized
            
            expectedNormMasses = [0.3; 0.7];
            expectedBinEdges = [0.0; 10.0; 20.0];
            expectedBinWidths = [10.0; 10.0];
            expectedPdfHeights = [0.03; 0.07];

            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges, expectedBinEdges, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths, expectedBinWidths, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(heights, expectedPdfHeights, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(sum(heights .* widths), 1.0, 'AbsTol', testCase.AbsTol);
        end
        
        function testThreeValuesUniformSpacing(testCase)
            % Test case 3: Three values, uniform spacing
            values = [10.0; 20.0; 30.0];
            probabilities = [1.0; 2.0; 1.0]; % Sum = 4
            
            expectedNormMasses = [0.25; 0.5; 0.25];
            expectedBinEdges = [5.0; 15.0; 25.0; 35.0];
            expectedBinWidths = [10.0; 10.0; 10.0];
            expectedPdfHeights = [0.025; 0.05; 0.025];

            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges, expectedBinEdges, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths, expectedBinWidths, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(heights, expectedPdfHeights, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(sum(heights .* widths), 1.0, 'AbsTol', testCase.AbsTol);
        end

        function testThreeValuesNonUniformSpacing(testCase)
            % Test case 4: Three values, non-uniform spacing
            values = [0.0; 2.0; 10.0];
            probabilities = [2.0; 2.0; 4.0]; % Sum = 8
            
            expectedNormMasses = [0.25; 0.25; 0.5];
            expectedBinEdges = [-1.0; 1.0; 6.0; 14.0];
            expectedBinWidths = [2.0; 5.0; 8.0];
            expectedPdfHeights = [0.125; 0.05; 0.0625];

            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges, expectedBinEdges, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths, expectedBinWidths, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(heights, expectedPdfHeights, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(sum(heights .* widths), 1.0, 'AbsTol', testCase.AbsTol);
        end
        
        function testProbabilitiesAreZero(testCase)
            % Test case 5: Probabilities sum to zero
            values = [10.0; 20.0; 30.0];
            probabilities = [0.0; 0.0; 0.0];
            
            expectedNormMasses = [0.0; 0.0; 0.0];
            expectedBinEdges = [5.0; 15.0; 25.0; 35.0]; 
            expectedBinWidths = [10.0; 10.0; 10.0];   
            expectedPdfHeights = [0.0; 0.0; 0.0];

            testCase.verifyWarning(@() vis.bayes.assess.pieceWiseConstantPDF(values, probabilities), ...
                                   'pieceWiseConstantPDF:ZeroTotalProbability');
            
            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges, expectedBinEdges, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths, expectedBinWidths, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(heights, expectedPdfHeights, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(sum(heights .* widths), 0.0, 'AbsTol', testCase.AbsTol);
        end

        function testExtremelyCloseValues(testCase)
            % Test case 6: Input values are distinct but extremely close
            funcTolerance = sqrt(eps); % Mirror the tolerance used inside the function
            delta = funcTolerance / 2.0; % Ensure delta is less than or equal to the function's tolerance
            
            values = [10.0; 10.0 + delta];
            probabilities = [1.0; 1.0]; % Sum = 2.0
            
            expectedNormMasses = [0.5; 0.5];
            expectedBinEdges = [values(1) - delta/2.0; values(1) + delta/2.0; values(2) + delta/2.0];
            expectedBinWidths = [delta; delta];
            expectedPdfHeights = [0.0; 0.0]; % Because widths <= funcTolerance
            expectedSumArea = 0.0;

            testCase.verifyWarning(@() vis.bayes.assess.pieceWiseConstantPDF(values, probabilities), ...
                                   'pieceWiseConstantPDF:ZeroWidthBin');

            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges, expectedBinEdges, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths, expectedBinWidths, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(heights, expectedPdfHeights, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(sum(heights .* widths), expectedSumArea, 'AbsTol', testCase.AbsTol);
        end
        
        function testLargeNumberOfPoints(testCase)
            % Test case 7: Larger number of points
            numPoints = 100;
            values = (1:numPoints)';
            probabilities = ones(numPoints, 1); 
            
            expectedNormMasses = ones(numPoints,1) / numPoints;
            
            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(values, probabilities);

            testCase.verifySize(normMasses, [numPoints, 1]);
            testCase.verifySize(widths, [numPoints, 1]);
            testCase.verifySize(edges, [numPoints+1, 1]);
            testCase.verifySize(heights, [numPoints, 1]);
            
            testCase.verifyEqual(normMasses, expectedNormMasses, 'AbsTol', testCase.AbsTol);
            testCase.verifyTrue(all(widths > 0));
            testCase.verifyEqual(sum(heights .* widths), 1.0, 'AbsTol', testCase.AbsTol);
            
            testCase.verifyEqual(edges(1), 0.5, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges(end), 100.5, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(edges(2), 1.5, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths(1), 1.0, 'AbsTol', testCase.AbsTol);
            testCase.verifyEqual(widths(2), 1.0, 'AbsTol', testCase.AbsTol);
        end

        % --- Tests for Input Validation and Error Handling ---
        function testValuesNotStrictlyIncreasingError(testCase)
            valuesFlat = [10.0; 10.0]; 
            probabilitiesFlat = [1.0; 1.0];
            testCase.verifyError(@() vis.bayes.assess.pieceWiseConstantPDF(valuesFlat, probabilitiesFlat), ...
                'validation:mustBeStrictlyIncreasing');

            valuesDecreasing = [10.0; 5.0]; 
            probabilitiesDec = [1.0; 1.0];
            testCase.verifyError(@() vis.bayes.assess.pieceWiseConstantPDF(valuesDecreasing, probabilitiesDec), ...
                'validation:mustBeStrictlyIncreasing');
        end

        function testValuesNotVectorError(testCase)
            valuesMatrix = [10.0 20.0; 30.0 40.0]; 
            probabilitiesMatrix = [1.0 1.0; 1.0 1.0]; 
             testCase.verifyError(@() vis.bayes.assess.pieceWiseConstantPDF(valuesMatrix, probabilitiesMatrix), ...
                'MATLAB:validation:IncompatibleSize');
        end

        function testProbabilitiesSizeMismatchError(testCase)
            values = [10.0; 20.0];
            probabilities = [1.0]; 
            testCase.verifyError(@() vis.bayes.assess.pieceWiseConstantPDF(values, probabilities), ...
                'validation:mustHaveSameSize');
        end
        
        function testNegativeProbabilitiesError(testCase)
            values = [10.0; 20.0];
            probabilities = [1.0; -1.0]; 
            testCase.verifyError(@() vis.bayes.assess.pieceWiseConstantPDF(values, probabilities), ...
                'MATLAB:validators:mustBeNonnegative');
        end
        
        function testEmptyInputsHandled(testCase)
            valuesEmpty = double.empty(0,1);
            probabilitiesEmpty = double.empty(0,1);
            
            testCase.verifyWarning(@() vis.bayes.assess.pieceWiseConstantPDF(valuesEmpty, probabilitiesEmpty), ...
                                   'MATLAB:pieceWiseConstantPDF:EmptyInput');
            
            [normMasses, widths, edges, heights] = vis.bayes.assess.pieceWiseConstantPDF(valuesEmpty, probabilitiesEmpty);
            testCase.verifyEmpty(normMasses);
            testCase.verifyEmpty(widths);
            testCase.verifyEmpty(edges);
            testCase.verifyEmpty(heights);
        end
    end
end
