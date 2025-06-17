function [IntensityChange, excludedCellsFig] = RelativeIntensityChange(MatchedCells)

%% Extract variables needed from the struct created in the matchCells function

matchedCentroids1 = MatchedCells.matchedCentroidCoordinates_0minus;
matchedRadii1 = MatchedCells.matchedRadii_0minus;
matchedCentroids2 = MatchedCells.matchedCentroidCoordinates_timepoint2; % Original Image2 coordinates
matchedRadii2 = MatchedCells.matchedRadii_timepoint2; % From original Image2
tform = MatchedCells.transformMatrix;
image1 = MatchedCells.GCaMP_Image1; % always 0minus (baseline)
image2 = MatchedCells.GCaMP_Image2; % 0plus in this case

%% Memory-Efficient Parallel Intensity Measurement
% Inner and outer parameters were optimized by systematically testing
% different combinations and validating with sensitivity analyses to
% balance data quality and quantity (only 0.19% intensity measurements
% invalid)
bgInnerFactor = 1.8;
bgOuterFactor = 2;

% Convert images to single precision to save memory
image1 = im2single(image1);
image2 = im2single(image2);

% Initialize arrays
% numMatched = sum(validMatches);
numMatched = size(matchedCentroids1,1);
intensity1 = zeros(numMatched, 1, 'single');
intensity2 = zeros(numMatched, 1, 'single');

% Precompute grid coordinates locally (do NOT broadcast full grids)
parfor i = 1:numMatched
    % --- Image 1 ---
    x1 = matchedCentroids1(i,1);
    y1 = matchedCentroids1(i,2);
    r1 = matchedRadii1(i);
    
    % Calculate local window coordinates (avoids full meshgrid)
    [rows, cols] = size(image1);
    [X_local1, Y_local1] = meshgrid(...
        max(1, floor(x1 - r1*bgOuterFactor)) : min(cols, ceil(x1 + r1*bgOuterFactor)), ...
        max(1, floor(y1 - r1*bgOuterFactor)) : min(rows, ceil(y1 + r1*bgOuterFactor)) ...
    );
    
    % Compute masks within the local window
    distSq1 = (X_local1 - x1).^2 + (Y_local1 - y1).^2;
    maskSoma1 = distSq1 <= r1^2;
    maskBG1 = distSq1 > (r1*bgInnerFactor)^2 & distSq1 <= (r1*bgOuterFactor)^2;
    
    % Extract local image patch
    localPatch1 = image1(...
        max(1, floor(y1 - r1*bgOuterFactor)) : min(rows, ceil(y1 + r1*bgOuterFactor)), ...
        max(1, floor(x1 - r1*bgOuterFactor)) : min(cols, ceil(x1 + r1*bgOuterFactor)) ...
    );
    
    intensity1(i) = median(localPatch1(maskSoma1)) - median(localPatch1(maskBG1));
    
    % --- Image 2 --- (Same approach)
    x2 = matchedCentroids2(i,1);
    y2 = matchedCentroids2(i,2);
    r2 = matchedRadii2(i);
    
    [rows2, cols2] = size(image2);
    [X_local2, Y_local2] = meshgrid(...
        max(1, floor(x2 - r2*bgOuterFactor)) : min(cols2, ceil(x2 + r2*bgOuterFactor)), ...
        max(1, floor(y2 - r2*bgOuterFactor)) : min(rows2, ceil(y2 + r2*bgOuterFactor)) ...
    );
    
    distSq2 = (X_local2 - x2).^2 + (Y_local2 - y2).^2;
    maskSoma2 = distSq2 <= r2^2;
    maskBG2 = distSq2 > (r2*bgInnerFactor)^2 & distSq2 <= (r2*bgOuterFactor)^2;
    
    localPatch2 = image2(...
        max(1, floor(y2 - r2*bgOuterFactor)) : min(rows2, ceil(y2 + r2*bgOuterFactor)), ...
        max(1, floor(x2 - r2*bgOuterFactor)) : min(cols2, ceil(x2 + r2*bgOuterFactor)) ...
    );
    
    intensity2(i) = median(localPatch2(maskSoma2)) - median(localPatch2(maskBG2));
    
end

% Calculate the change in calcium intensity for each cell
cell_FoldChange = intensity2./intensity1;
cell_log2FC = log2(cell_FoldChange);

% Find indices with significant imaginary components
significantImag = abs(imag(cell_log2FC)) > 1e-5;
problematicCells = find(significantImag);

cell_logFoldChange = cell_log2FC(~significantImag);
cell_logFoldChange(isinf(cell_logFoldChange)) = NaN;

% Plot the spatial distribution of excluded cells
excludedCellsFig = figure;
scatter(matchedCentroids1(:,1), matchedCentroids1(:,2), 'b', 'filled');
hold on;
scatter(matchedCentroids1(problematicCells,1), matchedCentroids1(problematicCells,2), 'r', 'filled');
legend('Valid', 'Excluded');
title('Spatial Distribution of Excluded Cells');


%% Save relevant variables into a structure
IntensityChange.cell_FoldChange = cell_FoldChange;
IntensityChange.cell_logFoldChange = cell_logFoldChange;
IntensityChange.cells_invalidLogs = problematicCells;
IntensityChange.mean_logFoldChange = mean(cell_logFoldChange,'omitnan');

% end

%% Check log2 calculations containing significant imaginary component
% % DETERMINE WHETHER TO PERFORM MEAN OR MEDIAN INTENSITY SUBTRACTION
% 
% % Load background pixel values array
% [bg1,bg1_loc] = uigetfile('*.mat','In 0minus folder - select bg.mat');
% [bg2,bg2_loc] = uigetfile('*.mat','In 0plus folder - select bg.mat');
% background1 = importdata([bg1_loc,bg1]);
% background2 = importdata([bg2_loc,bg2]);
% 
% % Find indices with significant imaginary components
% significantImag = abs(imag(cell_log2FC)) > 1e-5;
% problematicCells = find(significantImag);
% 
% % For each problematic cell, check the raw intensities and background regions:
% for i = 1:length(problematicCells)
%     idx = problematicCells(i);
% 
%     % Retrieve pre- and post-stretch intensities
%     preIntensity = intensity1(idx);
%     postIntensity = intensity2(idx);
% 
%     % Retrieve background regions
%     bgPre = background1(idx); % Cell array of background pixel values (pre-stretch)
%     bgPost = background2(idx); % Cell array of background pixel values (post-stretch)
% 
%     % Calculate mean/median background
%     bgMean_pre = mean(bgPre);
%     bgMedian_pre = median(bgPre);
%     bgMean_post = mean(bgPost);
%     bgMedian_post = median(bgPost);
% 
%     % Check if background subtraction caused negative intensities
%     if (preIntensity <= 0) || (postIntensity <= 0)
%         fprintf('Cell %d: Negative intensity detected (Pre: %.2f, Post: %.2f)\n', idx, preIntensity, postIntensity);
%     end
% end
% 
% % For each problematic cell, check the raw intensities and background regions:
% 
% % For a problematic cell (e.g., idx = 1)
% idx = problematicCells(1);
% 
% % Plot pre-stretch image with soma and background
% figure;
% imshow(image1);
% hold on;
% plot(matchedCentroids1(idx,1), matchedCentroids1(idx,2), 'ro', 'MarkerSize', 10);
% title(sprintf('Pre-Stretch: Cell %d (Intensity=%.2f)', idx, intensity1(idx)));
% 
% % Overlay background annulus
% theta = 0:0.1:2*pi;
% r1 = matchedRadii1(idx) * bgInnerFactor; % Inner radius
% r2 = matchedRadii1(idx) * bgOuterFactor; % Outer radius
% x_inner = r1 * cos(theta) + matchedCentroids1(idx,1);
% y_inner = r1 * sin(theta) + matchedCentroids1(idx,2);
% x_outer = r2 * cos(theta) + matchedCentroids1(idx,1);
% y_outer = r2 * sin(theta) + matchedCentroids1(idx,2);
% plot(x_inner, y_inner, 'g-', 'LineWidth', 1);
% plot(x_outer, y_outer, 'g-', 'LineWidth', 1);
% hold off;
% 
% % For each problematic cell, compare the mean and median background values:
% 
% % For a problematic cell
% idx = problematicCells(1);
% bgPre = background1(idx);
% bgPost = background2(idx);
% 
% fprintf('Cell %d:\n', idx);
% fprintf('Pre-Stretch Background: Mean=%.2f, Median=%.2f\n', mean(bgPre), median(bgPre));
% fprintf('Post-Stretch Background: Mean=%.2f, Median=%.2f\n', mean(bgPost), median(bgPost));
% 
% % Plot background pixel distribution
% figure;
% subplot(1,2,1);
% histogram(bgPre, 'BinWidth', .1);
% title('Pre-Stretch Background Pixels');
% subplot(1,2,2);
% histogram(bgPost, 'BinWidth', .1);
% title('Post-Stretch Background Pixels');
% dummy = 1; % dummy variable to stop code at to visualize the figures
end