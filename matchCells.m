function [MatchedCells, cellMatchesFig] = matchCells(folderpath, image1_folder, image2_folder, tform)

% Function to match the cells identified by NeuroCa between multiple
% timepoints for one sample (currently written for 4 timepoints)

%% Load data from NeuroCa calculations

% Load data from first timepoint
centroids1 = importdata([folderpath,image1_folder,'\center.mat']);
radii1 = importdata([folderpath,image1_folder,'\radii.mat']);
fdata1 = importdata([folderpath,image1_folder,'\fdata.mat']);

% Load data from the second timepoint being compared
centroids2 = importdata([folderpath,image2_folder,'\center.mat']);
radii2 = importdata([folderpath,image2_folder,'\radii.mat']);
fdata2 = importdata([folderpath,image2_folder,'\fdata.mat']);
if strcmp(image2_folder,'\1hr')
    image1 = im2double(imread([folderpath,image1_folder,'\beforeraw.tif']));
    image2 = im2double(imread([folderpath,image2_folder,'\firstimage.tif']));
elseif strcmp(image2_folder,'\24hr')
    image1 = im2double(imread([folderpath,image1_folder,'\beforeraw.tif']));
    image2 = im2double(imread([folderpath,image2_folder,'\firstimage.tif']));
elseif strcmp(image2_folder,'\0plus') && ~isempty(tform)
    image1 = im2double(imread([folderpath,image1_folder,'\beforeraw.tif']));
    image2 = im2double(imread([folderpath,image2_folder,'\afterraw.tif']));
end

%% Determine which frames from 0minus and 0plus to compare (immediately before and after stretch, excluding during transients)

if isempty(tform) && strcmp(image2_folder,'\0plus')
    
    % 0-minus
    fig1 = figure; clf;
    x_axis1 = fdata1(1,:);

    ax1 = subplot(5,1,1);
    plot(ax1, x_axis1,fdata1(3,:),'LineWidth',3);
    xlabel('Time (s)')
    title('SELECT ONE POINT RIGHT BEFORE STRETCH NOT WITHIN A TRANSIENT')
   
    ax2 = subplot(5,1,2);
    plot(ax2, x_axis1,fdata1(5,:),'LineWidth',3);
    xlabel('Time (s)')
    
    ax3 = subplot(5,1,3);
    plot(ax3, x_axis1,fdata1(10,:),'LineWidth',3);
    ylabel('Normalized Intensity')
    xlabel('Time (s)')
    
    ax4 = subplot(5,1,4);
    plot(ax4,x_axis1,fdata1(50,:),'LineWidth',3)
    xlabel('Tine (s)')

    ax5 = subplot(5,1,5);
    plot(ax5,x_axis1,fdata1(100,:),'LineWidth',3)
    xlabel('Time (s)')

    % Select desired x-axis value
    allAxes1 = [ax1 ax2 ax3 ax4 ax5];
    allY1 = {fdata1(3,:), fdata1(5,:), fdata1(10,:), fdata1(50,:), fdata1(100,:)};
    [xsel1, ~] = selectOnePoint(fig1, allAxes1, x_axis1, allY1);

    [~, frame_before] = min(abs(fdata1(1,:) - xsel1));

    % 0-plus
    fig2 = figure; clf;
    x_axis2 = fdata2(1,:);
    ax1 = subplot(5,1,1);
    plot(ax1, x_axis2,fdata2(3,:),'LineWidth',3);
    xlabel('Time (s)')
    title('SELECT ONE POINT RIGHT AFTER STRETCH NOT WITHIN A TRANSIENT')
   
    ax2 = subplot(5,1,2);
    plot(ax2, x_axis2,fdata2(5,:),'LineWidth',3);
    xlabel('Time (s)')
    
    ax3 = subplot(5,1,3);
    plot(ax3, x_axis2,fdata2(10,:),'LineWidth',3);
    ylabel('Normalized Intensity')
    xlabel('Time (s)')
    
    ax4 = subplot(5,1,4);
    plot(ax4,x_axis2,fdata2(50,:),'LineWidth',3)
    xlabel('Tine (s)')

    ax5 = subplot(5,1,5);
    plot(ax5,x_axis2,fdata2(100,:),'LineWidth',3)
    xlabel('Time (s)')
    
    % Select desired x-axis value
    allAxes2 = [ax1 ax2 ax3 ax4 ax5];
    allY2 = {fdata2(3,:), fdata2(5,:), fdata2(10,:), fdata2(50,:), fdata2(100,:)};
    [xsel2, ~] = selectOnePoint(fig2, allAxes2, x_axis2, allY2);

    [~, frame_after] = min(abs(fdata2(1,:) - xsel2));

       
    fprintf('\nSave "before" image at frame %2.0f \nSave "after" image + %2.0f frames from the starting 0plus frame previously determined \n\n',frame_before,frame_after-1)
    
    close all
    
    prompt = '*ENTER 1* when the two images ("beforeraw.tif" and "afterraw.tif") are saved in their respective before and after folders: ';
    continue_run = input(prompt);

    % Load original images
    image1 = im2double(imread([folderpath,image1_folder,'\beforeraw.tif']));
    image2 = im2double(imread([folderpath,image2_folder,'\afterraw.tif']));
end



%% Create transform matrices by manually matching point between the 0minus image and all other timepoints
if isempty(tform)
    %% --- Step 1: Manual Point Match ---
    % Get matching points from 0minus and 0plus images to account for any shifts between
    % before and after images
    numPoints = 6;  % Minimum 3 points required
    fprintf('Select %d points in Image 1 (left) first, then %d in Image 2 (right).\n', numPoints, numPoints);
    [points1, points2] = cpselect(image1, image2, 'Wait', true);

   
    % --- Step 2: Compute Transform & Validate ---
    tform = fitgeotrans(points2, points1, 'affine');

    % Calculate reprojection error on MANUAL POINTS (not centroids)
    points2_aligned = transformPointsForward(tform, points2);  % Transform manual points
    error_stretch = mean(sqrt(sum((points1 - points2_aligned).^2, 2)));
    fprintf('Mean alignment error (0minus-0plus): %.2f pixels\n', error_stretch);
    if error_stretch > 5
        warning('High alignment error between 0minus and 0plus - recheck manual points.');
    end
    
end

%% --- Step 3: Transform images and centroids ---
image2_aligned = imwarp(image2, tform, 'OutputView', imref2d(size(image1)));
centroids2_aligned = transformPointsForward(tform, centroids2);  % Applies to all centroids

%% --- Step 4: Match Somas ---
% maxDist = 10;
maxDist = 5;

% Comparing timepoints 1 and 2
% idx = sample 2 cell indices identified as the closest neighbor to each
% cell in sample 1
% dist = distance between the identified neighbors
[idx, dist] = knnsearch(centroids2_aligned, centroids1, 'K', 1);
validMatches = dist < maxDist; % only neighbors within the set max distance will be counted

% Extract matched pairs
matchedCentroids1 = centroids1(validMatches, :);
findIndices1 = 1:size(centroids1, 1);
matchedIndices(:,1) = findIndices1(validMatches)'; % cell numbers from 0minus

matchedIndices(:,2) = idx(validMatches); % Indices in ORIGINAL centroids2
matchedCentroids2 = centroids2(matchedIndices(:,2), :); % Original Image2 coordinates
matchedRadii1 = radii1(validMatches);
matchedRadii2 = radii2(matchedIndices(:,2)); % From original Image2

%% --- Step 5: Visualize Alignment ---
% Visualize alignment and matched somas
cellMatchesFig = figure;
imshowpair(image1, image2_aligned, 'falsecolor');
hold on;

% Plot Image1 centroids (before stretch)
plot(matchedCentroids1(:,1), matchedCentroids1(:,2), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 1.5);

% Plot Image2 centroids (after stretch, aligned)
% Use idx(validMatches) to get correct indices into centroids2_aligned
matchedIndices2 = idx(validMatches);
plot(centroids2_aligned(matchedIndices2,1), centroids2_aligned(matchedIndices2,2), 'g+', ...
    'MarkerSize', 10, 'LineWidth', 1.5);

% Draw lines between matched pairs
for i = 1:size(matchedCentroids1,1)
    plot([matchedCentroids1(i,1), centroids2_aligned(matchedIndices2(i),1)], ...
         [matchedCentroids1(i,2), centroids2_aligned(matchedIndices2(i),2)], ...
         'b-', 'LineWidth', 1);
end

hold off;
title(['Aligned Images and Soma Centroids (',num2str(size(matchedCentroids1,1)),' Cells Matched)']);
legend('0minus', [image2_folder(2:end),' (Aligned)'], 'Matches');

%% Save variables from matched pairs

MatchedCells.matchedCellNumbers = matchedIndices;
MatchedCells.matchedCentroidCoordinates_0minus = matchedCentroids1;
MatchedCells.matchedRadii_0minus = matchedRadii1;
MatchedCells.matchedCentroidCoordinates_timepoint2 = matchedCentroids2; % Original Image2 coordinates
MatchedCells.matchedRadii_timepoint2 = matchedRadii2; % From original Image2
MatchedCells.transformMatrix = tform;
MatchedCells.validMatches_image1 = validMatches; % logical whether each cell in sample 1 found a valid match
MatchedCells.GCaMP_Image1 = image1;
MatchedCells.GCaMP_Image2 = image2;

%% % --- Nested helper function ---
    function [xsel, idx] = selectOnePoint(figHandle, axesArray, xVals, yDataCells)
        % Initialize cursor line, dots, and labels
        for i = 1:numel(axesArray)
            ax = axesArray(i);
            hold(ax, 'on');
            % Cursor vertical line
            lines(i) = line(ax, [NaN NaN], ylim(ax), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
            % Red dot marker
            dots(i) = line(ax, NaN, NaN, 'Marker', 'o', 'Color', 'r', 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle', 'none');
            % Text label
            texts(i) = text(ax, NaN, NaN, '', 'Color', 'r', 'FontSize', 10, 'VerticalAlignment', 'bottom');
        end

        % Set mouse motion callback to update cursor line
        set(figHandle, 'WindowButtonMotionFcn', @(~,~) trackCursor(axesArray, lines));

        disp('Click on a point, then press Enter to confirm.');

        % Wait for user click
        [xsel, ~] = ginput(1);

        % Stop tracking cursor line
        set(figHandle, 'WindowButtonMotionFcn', '');

        % Find nearest index
        [~, idx] = min(abs(xVals - xsel));

        % Update cursor line, dots and labels on each axes
        for i = 1:numel(axesArray)
            ax = axesArray(i);
            yVals = yDataCells{i};
            set(lines(i), 'XData', [xVals(idx) xVals(idx)], 'YData', ylim(ax));
            set(dots(i), 'XData', xVals(idx), 'YData', yVals(idx));
            set(texts(i), 'Position', [xVals(idx), yVals(idx)], 'String', sprintf('%.2f s', xVals(idx)));
        end

        disp('Press Enter to continue.');
        waitforbuttonpress; % wait for Enter key
    end

end
%% Function to select the x-value using a vertical line cursor
function trackCursor(allAxes, allLines)
    cp = get(gcf, 'CurrentPoint');  % current point in figure units
    ax = gca;
    pt = get(ax, 'CurrentPoint');   % current point in data units

    xval = pt(1,1);  % extract x
    for i = 1:3
        if isvalid(allLines(i))
            ylimNow = ylim(allAxes(i));
            set(allLines(i), 'XData', [xval xval], 'YData', ylimNow);
        end
    end
end