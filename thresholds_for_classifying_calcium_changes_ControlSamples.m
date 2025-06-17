%% Use percentiles to define the threshold (Q-Q plot shows that the data is skewed)

% ------- READ BEFORE RUNNING THIS CODE -------
% After running baseline intensity change main code for all control/sham
% samples (depending on experiment), 
% *save IntensityChange.cell_logFoldChange as
% controlSamples_log2change.sample#
% *** save this new variable ('controlSamples_log2change') into the main
% directory containing all of your data

% Use control samples to define thresholds for classifying calcium
% fluorescent changes before and after stretch as "increase," "decrease,"
% and "no change"

log2FC_controls = [controlSamples_log2change.sample1; controlSamples_log2change.sample2; controlSamples_log2change.sample3; controlSamples_log2change.sample4; controlSamples_log2change.sample5; controlSamples_log2change.sample6; controlSamples_log2change.sample7; controlSamples_log2change.sample8; controlSamples_log2change.sample9];% controlSamples_log2change.sample10];
all_controls = log2FC_controls(~isnan(log2FC_controls));

% Compute percentile thresholds from controls
alpha = 2.5;
threshold_low = prctile(all_controls, alpha);
threshold_high = prctile(all_controls, 100 - alpha);

save('Control_Intensity_Change_Thresholds.mat', 'threshold_high','threshold_low','controlSamples_log2change')

%%
% Every mean of the log2(fold change) of the 4 control samples have
% been found to be <10% different of each other
% The sd are within 7% and 30% of each other
all_controls = figure;
hold on;
histogram(controlSamples_log2change.sample1, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample2, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample3, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample4, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample5, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample6, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample7, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample8, 'BinWidth', 0.05, 'Normalization', 'pdf');
histogram(controlSamples_log2change.sample9, 'BinWidth', 0.05, 'Normalization', 'pdf');
%histogram(controlSamples_log2change.sample10, 'BinWidth', 0.05, 'Normalization', 'pdf');

xlabel('log2(intensity2/intensity1)');
ylabel('Density');
title('Control Sample Distributions');
legend('Sample 1', 'Sample 2', 'Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9');%, 'Sample 10');

savefig(all_controls,'all_shams_log_foldChange.fig')
saveas(all_controls,'all_shams_log_foldChange.png')