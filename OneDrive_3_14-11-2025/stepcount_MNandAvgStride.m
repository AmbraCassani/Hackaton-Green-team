function [total_steps, steps_per_bout, cadence_spm, avg_stride_duration_table_results] = stepcount_MNandAvgStride(BA_ap, fs, bouts)
% This function computes total steps, step amount for each bout, cadence and average stride duration from AP trunk acceleration via unbiased autocorr.
% In the Moe-Nilssen (2003) paper, it is stated that the first positive lag peak of unbiased autocorrelation represents one step period, and second positive lag peak of unbiased autocorrelation represents one stride period.
% In the Moe-Nilssen (2003) paper, these equations are given: Number of steps M = N/n; cadence c = 60*fs/n.
% In the Weiss (2013) paper, step-frequency window used is 0.5–3 Hz.

detrended_BA = detrend(BA_ap(:));
%dur_min = max(1, ceil(fs/3)); dur_max = floor(fs/0.5);              % Max and min values of the step period range
[b,a] = butter(2, [0.5 3]/(fs/2), 'bandpass');
filtered_BA = filtfilt(b,a,detrended_BA);

num_bouts = numel(bouts);                                           % Number of valid bouts

% Initialization of variables
steps_per_bout = zeros(num_bouts,1);
cadence_spm = nan(num_bouts,1);
avg_stride_duration_per_bout = zeros(num_bouts,1);
total_steps = 0;


for b = 1:num_bouts
    start_bout = bouts{b}(1); 
    end_bout = bouts{b}(2);                             % Indices of the start and end of the bouts
    current_bout_signal = filtered_BA(start_bout:end_bout);
    N = numel(current_bout_signal);                                 % Number of samples in current bout

    Acorr = xcov(current_bout_signal,'unbiased');                                    % Unbiased autocorrelation applied to signal, which avoids attenuation across lags.
    Acorr = Acorr(ceil(length(Acorr)/2):end);                         % Taking only positive lags
    % lo = max(dur_min,1); 
    % hi = min(dur_max, numel(Apos));            % Lower and upper bounds of where we look for peaks in autocorrelation.
    % if hi <= lo
    %     continue; 
    % end

    [pks,locs] = findpeaks(Acorr);                            % Peaks from autocorrelation and their locations.
    % figure()
    % plot(1:length(Acorr),Acorr,locs,pks,'O')
    if isempty(pks), continue; end
    n = locs(1) - 1;                                           % Earliest strong peak's location in terms of samples.

    M = N / n;                                                      % Steps per bout M = N/n.
    avg_stride_duration_per_bout(b) = (locs(2)-1)/fs;   % Avg stride duration per bout [seconds]
    steps_per_bout(b) = max(0, round(M));                           % Rounded to integer steps
    cadence_spm(b)    = 60*fs / n;                                  % c = 60*fs/n as stated at the start.
    total_steps       = total_steps + steps_per_bout(b);
end
   avg_stride_duration_table_results = mean(avg_stride_duration_per_bout);
end
