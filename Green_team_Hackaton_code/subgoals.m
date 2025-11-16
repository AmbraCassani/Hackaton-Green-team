function [control, fallers] = subgoals(control, fallers, only_control, g, fs, win_sma, win_freqfilt, namelist, idx_co, idx_fl, threshold_control, threshold_fallers)


if only_control==1
    start_analysis = 1;
    end_analysis = 40;
elseif only_control == 0
        start_analysis = 41;
        end_analysis = length(namelist);
elseif only_control == 2
        start_analysis = 1;
        end_analysis = length(namelist);
end


for i=start_analysis:end_analysis
tic
    sgl = rdsamp(strcat('\long-term-movement-monitoring\',namelist{i}));
    sgl = sgl *g;
    [m,n] = size(sgl);
    t = (0:m-1)/fs;



    %%%%%%% The following steps are described in 11-appendix-3 %%%%%%%%%%


    [BA, sgl] = find_BA(sgl,[],fs,0);


    SMA = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);


    %%%%%%%%% PROCESSING AND BINARY SIGNAL OF FREQUENCY FILTER %%%%%%%%%%%%

    binary_signal_freqfilt = compute_frequency_activity_mask(sgl,fs,win_freqfilt);


    %%%%%%%%%%%%%%%%%%%%%%% BINARY SIGNAL SMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    binary_signal_SMA = sma_thresholding(sgl,threshold_control,threshold_fallers,SMA,namelist,i);


    %%%%%%%%%%%%%%%%%%%%%%%%% BINARY SIGNAL MERGING %%%%%%%%%%%%%%%%%%%%%%%
    binary_signal = double(logical(binary_signal_SMA) | logical(binary_signal_freqfilt));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    actv_duration_samples = extract_walking_bouts(binary_signal,fs);
  

    % Used to compute average stride duration and total number of steps for
    % each subject
    [total_steps(i), steps_per_bout, cadence_spm, avgstride] = stepcount_MNandAvgStride(BA(:,3),fs,actv_duration_samples);


    % Transform the boundary indices of activity into effective time period
    actv_duration_time = samples_to_time(actv_duration_samples,fs);

    total_walking_time = sum(actv_duration_time); % [seconds]
    total_walking_percent = total_walking_time*100/t(end);

    disp(i)



    [control, fallers, idx_co, idx_fl] = memorization_subgroup(control, fallers, namelist, total_walking_percent, total_steps(i), avgstride, i, idx_co, idx_fl);

    toc
end


% Results calculation
if only_control == 1 || only_control == 2
    control.meanwp = mean(control.walking_percent);
    control.stdwp = std(control.walking_percent);
    control.medianwp = median(control.walking_percent);
    idx_outl_wp = findOutliers(control.walking_percent);
    control.meanwpNOoutl = mean(control.walking_percent(idx_outl_wp));
    control.stdwpNOoutl = std(control.walking_percent(idx_outl_wp));

    control.meants = mean(control.total_steps);
    control.stdts = std(control.total_steps);
    control.mediants = median(control.total_steps);
    idx_outl_ts = findOutliers(control.total_steps);
    control.meantsNOoutl = mean(control.total_steps(idx_outl_ts));
    control.stdtsNOoutl = std(control.total_steps(idx_outl_ts));


    control.meanavgstride = mean(control.avgstride);
    control.stdavgstride = std(control.avgstride);
    control.medianavgstride = median(control.avgstride);
    idx_outl_avgstride = findOutliers(control.avgstride);
    control.meanavgstrideNOoutl = mean(control.avgstride(idx_outl_avgstride));
    control.stdavgstrideNOoutl = std(control.avgstride(idx_outl_avgstride));
end


if ~only_control || only_control == 2
    fallers.meanwp = mean(fallers.walking_percent);
    fallers.stdwp = std(fallers.walking_percent);
    fallers.medianwp = median(fallers.walking_percent);
    idx_outl_wp = findOutliers(fallers.walking_percent);
    fallers.meanwpNOoutl = mean(fallers.walking_percent(idx_outl_wp));
    fallers.stdwpNOoutl = std(fallers.walking_percent(idx_outl_wp));


    fallers.meants = mean(fallers.total_steps);
    fallers.stdts = std(fallers.total_steps);
    fallers.mediants = median(fallers.total_steps);
    idx_outl_ts = findOutliers(fallers.total_steps);
    fallers.meantsNOoutl = mean(fallers.total_steps(idx_outl_ts));
    fallers.stdtsNOoutl = std(fallers.total_steps(idx_outl_ts));


    fallers.meanavgstride = mean(fallers.avgstride);
    fallers.stdavgstride = std(fallers.avgstride);
    fallers.medianavgstride = median(fallers.avgstride);
    idx_outl_avgstride = findOutliers(fallers.avgstride);
    fallers.meanavgstrideNOoutl = mean(fallers.avgstride(idx_outl_avgstride));
    fallers.stdavgstrideNOoutl = std(fallers.avgstride(idx_outl_avgstride));
end






% Results visualization
if only_control == 1 || only_control == 2

    fprintf('\n===== RESULTS NON FALLERS =====\n');
    fprintf('The mean value of the total walking percentage is %.2f %%\n', control.meanwpNOoutl);
    fprintf('The standard deviation of the total walking percentage is %.2f \n', control.stdwpNOoutl);
    fprintf('The mean value of the total steps is %.2f %%\n', control.meantsNOoutl);
    fprintf('The standard deviation of the total steps is %.2f \n', control.stdtsNOoutl);

end

if ~only_control || only_control == 2
    fprintf('\n===== RESULTS FALLERS =====\n');
    fprintf('The mean value of the total walking percentage is %.2f %%\n', fallers.meanwpNOoutl);
    fprintf('The standard deviation of the total walking percentage is %.2f \n', fallers.stdwpNOoutl);
    fprintf('The mean value of the total steps is %.2f %%\n', fallers.meantsNOoutl);
    fprintf('The standard deviation of the total steps is %.2f \n', fallers.stdtsNOoutl);
end


end