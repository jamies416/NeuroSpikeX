function [Mean_Ca_change_fold, cell_Ca_change_fold, mean_raw_intensity, std_raw_intensity, mean_intensity_fig] = calciumRawIntensityChange(result, frame_rate)


% Loop through each time frame and average the raw intensity values of all the cells
for i = 1:size(result,2)
    
    mri = mean(result(2:end, i),'all','omitnan');
    mean_raw_intensity(i) = mri;

    ss1 = std(result(2:end,i),[],'all','omitnan');
    std_raw_intensity(i) = ss1; 
end

% Specifically compare the change in background intensity
f = figure;
x_axis = (1:size(result,2))./frame_rate;
plot(x_axis,result(2,:),'LineWidth',3);
ylabel('Raw Intensity (a.u.)')
xlabel('Time (s)')
xlim([-inf size(result,2)/frame_rate])
title('DRAG A RECTANGLE FROM RIGHT BEFORE THE STRETCH TO RIGHT AFTER')
[~,xdata,~] = selectdata('selectionmode','rect');
frame_before = xdata(1,1)*frame_rate;
frame_after = xdata(end,1)*frame_rate;

frames_5s = round(frame_rate*5); % # of frames in 5-seconds
background_before = result(2:end, (frame_before-frames_5s):frame_before); % background raw intensity of each cell before stretch
background_after = result(2:end, frame_after:(frame_after+frames_5s)); % background raw intensity of each cell after stretch

avg_background_before = mean(background_before,2); % average background before stretch for each cell
avg_background_after = mean(background_after,2); % average background after stretch for each cell

% Calculate the change in calcium intensity for each cell
cell_Ca_change_fold = avg_background_after./avg_background_before;

% Calculate average calcium intensity change for the entire network
Mean_Ca_change_fold = mean(cell_Ca_change_fold);

% close

%% PLOT MEAN RAW INTENSITY WITH SD

frame_num = 1:size(result,2);
time = frame_num./frame_rate;

mean_intensity_fig = figure;
shadedErrorBar(time, mean_raw_intensity, std_raw_intensity,{'-','Color',[0.00,0.45,0.74]},1)
xlabel('Time (s)')
ylabel('Raw Intensity (a.u.)')
xlim([0,max(time)]);
set(gca, 'FontName', 'Verdana')
title('Mean Raw Intensity Before and After Stretch')

%% SUMMARY OF "MIDDLE" DATA
summary = msgbox(sprintf('Calcium Raw Intensity Change (a.u. fold) = %2.3g',Mean_Ca_change_fold),'DATA DURING STRETCH SUMMARY', 'non-modal');

set(summary, 'position', [500 440 350 75]); %makes box bigger
ah = get( summary, 'CurrentAxes' );
ch = get( ah, 'Children' );
set( ch, 'FontSize', 14 ); %makes text bigger

end