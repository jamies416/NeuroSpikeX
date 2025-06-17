%%
clearvars
close all
clc

%%
folderpath = uigetdir(pwd,'*Select folder of entire sample dataset to analyze*');

%% Find cell matches between 0-minus and 0-plus/1-hr/24-hr

% If transforms have already been calculated, the last input into the
% following functions should be the transformation matrix
[MatchedCells_0minus_0plus, cellMatchesFig_0minus_0plus] = matchCells(folderpath,'\0minus','\0plus',[]);
% [MatchedCells_0minus_0plus, cellMatchesFig_0minus_0plus] = matchCells(folderpath,'\0minus','\0plus',MatchedCells_0minus_0plus.transformMatrix);

% [MatchedCells_0minus_1hr, cellMatchesFig_0minus_1hr] = matchCells(folderpath,'\0minus','\1hr',[]);
% [MatchedCells_0minus_24hrs, cellMatchesFig_0minus_24hr] = matchCells(folderpath,'\0minus','\24hr',[]);

% Save variables and figures
% save('MatchedCells.mat','MatchedCells_0minus_0plus','MatchedCells_0minus_1hr','MatchedCells_0minus_24hrs')
save('MatchedCells.mat','MatchedCells_0minus_0plus')
savefig(cellMatchesFig_0minus_0plus,'matchedCellsFig_0minus_0plus.fig')
saveas(cellMatchesFig_0minus_0plus,'matchedCellsFig_0minus_0plus.png')
% savefig(cellMatchesFig_0minus_1hr,'matchedCellsFig_0minus_1hr.fig')
% saveas(cellMatchesFig_0minus_1hr,'matchedCellsFig_0minus_1hr.png')
% savefig(cellMatchesFig_0minus_24hr,'matchedCellsFig_0minus_24hr.fig')
% saveas(cellMatchesFig_0minus_24hr,'matchedCellsFig_0minus_24hr.png')

%% Calculate Relative Fluorescent Intensity Change Immediately After Stretch

% Calculate the relative intensity change and log2(fold change) and plot
% matched cells excluded from this calculation
close all

[IntensityChange, excludedCellsFig] = RelativeIntensityChange(MatchedCells_0minus_0plus);

% After running all control/sham samples and the
% "thresholds_for_classifying_calcium_changes_ControlSamples.m" uncomment
% the following line
% Uncomment line below for strain experiments:
%load('\\franck-syn.me.wisc.edu\francklabdata\Individual\Jamie\Paper 1 Experiments\Strain Threshold - GCaMP6s\Sham\Control_Intensity_Change_Thresholds.mat')

% Uncomment line below for glutamate experiments:
load('\\franck-syn.me.wisc.edu\francklabdata\Individual\Jamie\Paper 1 Experiments\Glutamate - GCaMP6s\Control\Control_Intensity_Change_Thresholds.mat')

% Plot distributions of calcium baseline change
normalized_overallCaChangeHist = NormalizedCalciumDataHistogram(IntensityChange.cell_logFoldChange,0,0,0,1,'log_2(Intensity Fold)','a.u.',round(std(IntensityChange.cell_logFoldChange,'omitnan')/6,2),'Normal');
overallCaChangeHist = CalciumDataHistogram(IntensityChange.cell_logFoldChange,0,0,0,1,'log_2(Intensity Fold)','a.u.',round(std(IntensityChange.cell_logFoldChange,'omitnan')/6,2),'Normal');

save('IntensityChange_stretch.mat','IntensityChange')
savefig(excludedCellsFig,'excluded_cells_beforeafter.fig')
saveas(excludedCellsFig,'excluded_cells_beforeafter.png')

savefig(normalized_overallCaChangeHist,'normalized_overallCaChangeHistogram.fig')
saveas(normalized_overallCaChangeHist,'normalized_overallCaChangeHistogram.png')

savefig(overallCaChangeHist,'corrected_overallCaChangeHistogram.fig')
saveas(overallCaChangeHist,'corrected_overallCaChangeHistogram.png')



%% Visualize Spatial Distribution of log2(fold change)
% Change threshold determined by negative control data from 4 samples
matchedCentroids = MatchedCells_0minus_0plus.matchedCentroidCoordinates_0minus;

% Classify experimental cells
isIncrease = IntensityChange.cell_logFoldChange > threshold_high;
IncreaseColor = '#1F77B4';
isDecrease = IntensityChange.cell_logFoldChange < threshold_low;
DecreaseColor = '#D62728';
isNoChange = ~isIncrease & ~isDecrease;
sameColor = '#FFFF00';

log2foldSpatialDistribution = figure;
imshow(MatchedCells_0minus_0plus.GCaMP_Image1)
hold on
scatter(matchedCentroids(isIncrease,1), matchedCentroids(isIncrease,2), 'o', ...
    'MarkerFaceColor',IncreaseColor,'MarkerEdgeColor',IncreaseColor,'MarkerFaceAlpha',0.6,'SizeData',75, 'LineWidth', 1.5);
scatter(matchedCentroids(isNoChange,1), matchedCentroids(isNoChange,2), 'o', ...
    'MarkerFaceColor',sameColor,'MarkerEdgeColor',sameColor,'MarkerFaceAlpha',0.5, 'SizeData',75,'LineWidth', 1.1);
scatter(matchedCentroids(isDecrease,1), matchedCentroids(isDecrease,2), 'o', ...
    'MarkerFaceColor',DecreaseColor,'MarkerEdgeColor',DecreaseColor,'MarkerFaceAlpha',0.5, 'SizeData',75,'LineWidth', 1.5);

title('Baseline Calcium Fluorescent Change Due to Stretch')
lgd = legend(['Increase (',num2str(sum(isIncrease)),' cells)'], ['No Change (',num2str(sum(isNoChange)),' cells)'], ['Decrease (',num2str(sum(isDecrease)),' cells)']);
fontsize(lgd,16,'points')

savefig(log2foldSpatialDistribution,'overallCaChange_spatialDistribution.fig')
saveas(log2foldSpatialDistribution,'overallCaChange_spatialDistribution.png')

fprintf('\nMean log2(fold change) = %4.4f \n# Cells Increased = %d \n# Cells No Change = %d \n# Cells Decreased = %d\n',IntensityChange.mean_logFoldChange,sum(isIncrease),sum(isNoChange),sum(isDecrease))
close all