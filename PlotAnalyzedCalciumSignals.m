function [testCalciumSignal_figure] = PlotAnalyzedCalciumSignals(cell_num, smooth_signal, signal_info_struct, calcium_transients_struct, rate_constants_struct, prm)
%% Extract necessary variables from input structures

fdata = signal_info_struct.fdata; % all cells
active_cells = calcium_transients_struct.active_cells; % contains the cell number of all active cells
true_pks = calcium_transients_struct.true_peaks; % only active cells
pk_locs = calcium_transients_struct.peak_locs; % only active cells
min_locs = calcium_transients_struct.decay_min_locs; % only active cells
smooth_diff = calcium_transients_struct.smooth_diff; % only active cells
time_axis = signal_info_struct.time;
% excluded_cst = rate_constants_struct.excluded_rate_cst;

%%
last = size(smooth_signal,2); % final frame to plot of shifted data

% Return an error if the cell entered is not an active cell

if isempty(find(cell_num == active_cells(:,1), 1))
    testCalciumSignal_figure = figure;

    plot(time_axis(1,1:last), fdata(cell_num,1:last))
    hold on
    plot(time_axis(1,1:last), smooth_signal(cell_num,1:last), 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('dF/F')
    xlim([0 120])
    title(['Cell #',num2str(cell_num),' Classified Inactive: choose a different cell number'])
    hold off
    return 
end

active_cell_index = find(cell_num == active_cells(:,1), 1);
% prm = 0.7;

last2 = size(smooth_diff,2); % final frame to plot of shifted derivative

% find the values of the min_locs (on the smooth_signal)
total_num_min = size(min_locs(active_cell_index,:), 2 );
min_values = zeros(1,total_num_min);
excluded_mins_loc = zeros(1,total_num_min);
excluded_mins_value = zeros(1,total_num_min);

for i = 1 : total_num_min
    if isnan(min_locs(active_cell_index,i))
        continue
    end
    min_values(1,i) = smooth_signal(cell_num, min_locs(active_cell_index,i));

    % if isnan(excluded_cst(active_cell_index,i))
    %     continue
    % end
    % excluded_mins_loc(1,i) = min_locs(active_cell_index,i);
    % excluded_mins_value(1,i) = min_values(1,i);
end
% excluded_mins_loc(excluded_mins_loc==0) = NaN;
% excluded_mins_value(excluded_mins_value==0) = NaN;

% Plot
testCalciumSignal_figure = figure;
tiledlayout(3,1)

nexttile
plot(time_axis(1,1:last), fdata(cell_num,1:last))
hold on
plot(time_axis(1,1:last), smooth_signal(cell_num,1:last),'LineWidth',2)
xlim([0 120])
xlabel('Time (s)')
ylabel('dF/F')
title(['Filtered Signal Cell #', num2str(cell_num)])
hold off

% nexttile
% plot(time_axis(1,1:last-1),diff_avgsignal)
% title('Differentiate of Signal')
% hold on
% plot(time_axis(1,1:last3),smooth_diff)

nexttile
plot(smooth_diff(active_cell_index,1:last2))
hold on
plot(smooth_signal(cell_num,1:last2),'LineWidth',2)
legend('Smoothed Derivative','Smoothed Signal')
xlabel('Frame #')
title('Exponential moving average of signal derivative)')
hold off

nexttile
plot(smooth_signal(cell_num,:),'b-','LineWidth',2)
hold on
plot(pk_locs(active_cell_index,:),true_pks(active_cell_index,:),'g.','MarkerSize',30)
plot(min_locs(active_cell_index,1:total_num_min),min_values(1,:),'r.','MarkerSize',30)
% if ~isempty(find(~isnan(excluded_cst(active_cell_index,:))==1,1))
%     plot(excluded_mins_loc(1,1:total_num_min),excluded_mins_value(1,:),'xk','MarkerSize',20,'LineWidth',2)
%     legend('Filtered Signal', 'Transient Peak', 'Transient End','Excluded Decays')
% else
    legend('Filtered Signal', 'Transient Peak', 'Transient End')
% end
xlim([0 size(smooth_signal,2)])
title(['Calcium transient peaks (minimum prominence ',num2str(prm),')'])
xlabel('Frame #')
ylabel('dF/F')
fontsize (14,'points');
hold off

%% Find decay constant from outputs of the above 2 sections

% Use the transient peaks and ends to create a decay window for each
% calcium spike
% Use the window on the un-filtered df/f signal to fit an exponential

% Fit exponential. Get the decay constant for every transient in every cell
k = 1;
i = 1; %just for the sake of testing 1 cell
figure()
tcl = tiledlayout(ceil(sum(~isnan(true_pks(active_cell_index,:)))/4),4);
% tcl = tiledlayout(1,2);
title(tcl,['Fitted Calcium Spike Decays from Cell #',num2str(cell_num)])
for j = 1:size(true_pks,2) % for loop runs through each transient
    % if isnan(true_pks(active_cell_index,j)) || ~isnan(excluded_cst(active_cell_index,j))
    %     continue
    % end
    clear obj_fun sol x y y_fit coeffs
    peak_frame = pk_locs(active_cell_index,j);
    end_frame = min_locs(active_cell_index,j);

    if isnan(peak_frame)
        continue
    end
    
    % generate data window
    x = time_axis(1,1:(end_frame-peak_frame+1)); % shifting all exponentials to begin at time 0 
    y = fdata(cell_num,peak_frame:end_frame) - fdata(cell_num,end_frame);  % decay signal (starting at peak and shifted to end at y=0)

    % exponential fit method 
    [y_fit,gof] = fit(x', y', 'exp1');
    coeffs = coeffvalues(y_fit);
    R(i,1) = cell_num; % first column is the cell number
    R(i,k+1) = -coeffs(2); % The decay rate constants (fitted) for every signal in each cell
    % Note: larger rate constant means a faster decay
    k = k+1;

    % Plot fits on each spike decay
    nexttile
    plot(y_fit,'r-', x, y, 'b-')
    grid on;
    % hold on;
    % plot(x, y_fit, 'r-', 'LineWidth', 2);
    % grid on;
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('dF/F', 'FontSize', 14);
    legendHandle = legend('Signal', 'Exp. Fit', 'Location', 'northeast');
    legendHandle.FontSize = 12;
    % Place formula text roughly in the middle of the plot.
    formulaString = sprintf('Y = %.3f * exp(%.7f * X)', coeffs(1), coeffs(2));
    xl = xlim;
    yl = ylim;
    xt = xl(1) + abs(xl(2)-xl(1)) * 0.325;
    yt = yl(1) + abs(yl(2)-yl(1)) * 0.59;
    text(xt, yt, formulaString, 'FontSize', 12, 'FontWeight', 'bold');

end

% figure
% plot(R(1,2:end),'.','MarkerSize',35)
% grid on
% xlabel('Calcium Transient Number')
% ylabel('Decay Rate Constant (s^{-1})')
% title(['Cell #',num2str(cell_num)])

end