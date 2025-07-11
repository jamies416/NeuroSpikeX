function [signal_information] = ProcessNeuroCa(fdata_original, desired_duration_seconds)

%% Function to get rid of any Nan of Inf values in the fdata that is sometimes an artifact of the NeuroCa calculations.
% Also, pull out all relevant data needed for subsequent analysis
%
% INPUTS:
% fdata_original = the fdata matrix generated by NeuroCa 
%                  (first row = time axis, other rows = cells, columns = frames)
% desired_duration_seconds = the target duration for analysis in seconds
%                            (from settings saved in Preprocessing Settings App)
%
% OUTPUTS:
% time_axis = isolated time axis
% fdata = matrix of just cell signals, checked for any NaN or inf artifacts
%         generated by NeuroCa
% frame_rate = image acquisition frame rate for each sample
% f_end = final frame # to ensure equal signal lengths of 120 sec

%%
    % Pull out essential information from NeuroCa's fdata output (fdata_original)
    time_axis = fdata_original(1,:); % time axis created by NeuroCa
    fdata_original(1,:) = []; %delete first row, now fdata contains only cell signals
    fdata = fdata_original;

    num_cells_i = size(fdata,1); % # of cells

    % Get rid of NaN or inf rows in fdata to clean up data
    for i = 1:num_cells_i
        if isnan(fdata_original(i,:))
            fdata(i,:) = NaN;
        elseif isinf(fdata_original(i,:))
            fdata(i,:) = NaN;
        end
    end
    
    fdata(isnan(fdata))=0;
    
    % Delete the cells that have no signal
    gone = find(~sum(fdata,2)); % locates the rows of all 0's
    fdata(gone,:) = [];

    if isempty(fdata) || isempty(time_axis) || time_axis(end) <= 0 % Check for valid data
        warning('ProcessNeuroCa: fdata or time_axis is empty or invalid after initial processing. Cannot calculate frame rate or f_end.');
        signal_information.time = time_axis; % Return what we have
        signal_information.fdata = fdata;
        signal_information.fps = NaN;
        signal_information.analysis_window_frames = 0; % New field name
        return;
    end

    % Create analysis time window (specific to user's needs identified in
    % the Preprocessing Settings App)
    frame_rate = (size(fdata,2)-1)/(time_axis(end) - time_axis(1)); % frames per second
    
    
    if frame_rate <=0 || isnan(frame_rate) || isinf(frame_rate)
        warning('ProcessNeuroCa: Calculated frame rate is invalid (%.2f fps). Check time axis and fdata dimensions. Setting analysis_window_frames to full length.', frame_rate);
        frame_rate = NaN; % Mark as invalid
        analysis_window_frames_target = size(fdata,2); % Default to full length if fps is bad
    else
        analysis_window_frames_target = round(desired_duration_seconds * frame_rate);
    end
    
    % Determine the actual number of frames to use for the analysis window
    % It's the minimum of the target duration (in frames) and the actual signal length
    actual_analysis_window_frames = min(analysis_window_frames_target, size(fdata,2));
    
    if actual_analysis_window_frames < 1 % Ensure it's at least 1 if there's data
        actual_analysis_window_frames = size(fdata,2); % Fallback to full length if calculation is odd
        if actual_analysis_window_frames < 1, actual_analysis_window_frames = 0; end
    end

    %% SAVE VARIABLES INTO A STRUCTURE

    signal_information.time = time_axis;
    signal_information.fdata = fdata; % full original fdata length (not truncated)
    signal_information.fps = frame_rate;
    signal_information.analysis_window_frames = actual_analysis_window_frames; % frame at desired analysis time length
    
end