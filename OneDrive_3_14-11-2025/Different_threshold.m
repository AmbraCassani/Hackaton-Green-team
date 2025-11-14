clear all
close all
clc
tic
%% Uncomment and run this section the first time you lunch the program to install the wfdb toolbox
% [old_path]=which('rdsamp'); if(~isempty(old_path)) rmpath(old_path(1:end-8)); end
% wfdb_url='https://physionet.org/physiotools/matlab/wfdb-app-matlab/wfdb-app-toolbox-0-10-0.zip';
% [filestr,status] = urlwrite(wfdb_url,'wfdb-app-toolbox-0-10-0.zip');
% unzip('wfdb-app-toolbox-0-10-0.zip');

%% Preprocess

mcodepath = cd;
mcodepath = strcat(cd,'\mcode');
addpath(mcodepath); %used to not having to manually-re add the folder each run

%% Data reading and processing

direc = dir('long-term-movement-monitoring/LabWalks/');
direc(1:2) = [];
direc(2:2:end) = [];


namelist = {};
for i=1:length(direc)
    namelist{i,1} = direc(i).name(1:10); % 1:10 gives only the name without the format, usefull for the main for loop
    namevalue(i) = str2num(direc(i).name(4:5));
end



% missingDataidx = find(diff(namevalue)>1) +1;
% for i=1:length(missingDataidx)
%     namelist = { namelist{1:missingDataidx(i)-1} {'MissingData'} namelist{missingDataidx(i):end} };
% end
fs = 100;
win_sma = 1*fs;
D_moe = 2*25; %length of the long hallway
h_zejl=ones(length(namelist),1);
alpha = 0.5; % relation between leg length and height

idx_co = 1;
idx_fl = 1;

for i= 1:length(namelist)
    %sgl = rdsamp('ltmm/CO001'); %use this to download it from the remote server

    sgl = rdsamp(strcat('long-term-movement-monitoring/LabWalks/',namelist{i}));
    g = 9.81; %TILT CORRECTION METHOD
    sgl(:,1:3) = sgl(:,1:3) * g; %convert to m/s^2
    [m,n] = size(sgl);

    t = (0:m-1)/fs;
    %labelplot = ["Vertical acceleration [g]"  "Mediolateral acceleration [g]" "AnteroPosterior acceleration [g]" "Yaw acceleration [degrees/s]" ...
    %"Pitch acceleration [degrees/s]" "Roll acceleration [degrees/s]"];
    %tic
    % figure()
    % for i = 1:6
    %
    %     subplot(3,2,i)
    %     plot(t,sgl(:,i))
    %     xlabel('time [sec]')
    %     ylabel(labelplot(i))
    %     xlim([0,t(end)])
    %
    % end
    % toc

    %% Filtering and preprocessing
    % sgl = sgl-mean(sgl);
    % sgl = detrend(sgl);
    % sgl(:,1) = detrend(sgl(:,1));
    % [b, a] = butter(4,30/fs);
    % sgl(:,1:3) = filtfilt(b,a,sgl(:,1:3));
    % figure(1)
    % freqz(b,a,[],fs)
    %[~,~,sgl(:,1)] = algo_Moe_Nilssen(sgl(:,3),sgl(:,2),sgl(:,1),'tiltAndNoG'); %TILT CORRECTION METHOD
    %% Correlation and peak detection

    %sgl(:,1) = sgl(:,1)-mean(sgl(:,1));
    correlation = xcov(sgl(:,3),'unbiased');
    maxpeak = max(correlation);
    correl_scaled = correlation/maxpeak;
    % figure()
    % plot(1:length(correlation),correlation)
    % figure()
    % plot(1:length(correlation),correl_scaled)
    [p_amp,p_locs] = findpeaks(correl_scaled,'MinPeakProminence',0.2,'MinPeakHeight',0.1); %to avoid some weird results for cadence

    %% Adaptive threshold peakfinder
    window = round(0.5 * fs); % 0.5 s window
    mov_avg = movmean(correl_scaled, window);
    mov_std = movstd(correl_scaled, window);

    % Step 2: Adaptive threshold based on statistical variation
    k = 1.5; % Scaling constant
    adaptive_thresh = mov_avg +  mov_std;

    [p_amp,p_locs] = findpeaks(correl_scaled);


    p_locs_threshed = p_locs(p_amp>adaptive_thresh(p_locs));
    p_amp_threshed = correl_scaled(p_locs_threshed);


    % figure(2)
    % subplot(2,1,1)
    % plot(t,sgl(:,1))
    % xlim([0 t(end)])
    % title(namelist{i}(1:5))
    % subplot(2,1,2)
    % plot(1:length(correlation),correl_scaled,p_locs_threshed,p_amp_threshed,'O',1:length(correlation),adaptive_thresh,'LineWidth',1)
    % xlim([0 length(correlation)])
    % pause;

    p_locs = p_locs_threshed;


    %% Procedural plot
    % figure(2)
    % subplot(2,1,1)
    % plot(t,sgl(:,1))
    % xlim([0 t(end)])
    % title(namelist{i}(1:5))
    % subplot(2,1,2)
    % plot(1:length(correlation),correl_scaled,p_locs,p_amp,'O')
    % xlim([0 length(correlation)])
    % pause;

    %% Turning detection to improve SMA computation

    %setting thresholds: values reported in the article
    thresholdDegrees=45;%only turns with degrees equal or over this threshold are identified
    thresholdLPfiltering=1.5;%Hz
    thresholdPeakDS=15;%degrees/second
    thresholdCrossingDS=5;%degrees/second
    thresholdIntraTurnS=0.05;%seconds
    thresholdTurnDurationLowS=0.5;
    thresholdTurnDurationHighS=10;

    %[Report, myreportallturns, idx_turns] = myfindTurns(sgl(:,4),fs,15,5,0.05,0.54,10,90);
    [idx_turns] = findTurns_prof(sgl(:,4), fs, thresholdDegrees, ...
        thresholdLPfiltering, thresholdPeakDS,thresholdCrossingDS,thresholdIntraTurnS, ...
        thresholdTurnDurationLowS, thresholdTurnDurationHighS );
    %%%% END OF TURNING DETECTION



    n_peaks = length(p_locs);
    center = ceil(n_peaks/2);
    half_p_locs = p_locs(center:end);

    dom_period = diff(half_p_locs(1:2));
    % if dom_period> (60*fs)/180 %plausible value for n based on the cadence (180),this is done to exclude some abnormal results due to errors in correlation peak detection
    %     n_moe = dom_period; %samples/dominant_period
    % else
    %     n_moe = half_p_locs(3)-half_p_locs(1);
    % end
    n_moe = dom_period;

    S_moe = t(end); % [seconds] time to walk D meters
    N_moe = S_moe*fs; %sample/D meters => length of our signal
    M_moe = N_moe/n_moe; %steps/Dmeters



    %%%%%%%%%%%% Zejlstra %%%%%%%%%%%%%%%%%%%%%%
    sgl_zejl = sgl(:,1)-mean(sgl(:,1));
    sgl_zejl(idx_turns) = [];
    t_zejl = (0:length(sgl_zejl)-1)/fs;
    [b, a] = butter(4,20/(fs/2),"low");
    sgl_zejl = filtfilt(b,a,sgl_zejl);
    velocity = cumtrapz(t_zejl,sgl_zejl(:,1));
    % position = cumtrapz(t_zejl,velocity);
    % [b, a] = butter(4,0.5/(fs/2),"high");
    % pos_filtered = filtfilt(b,a,position);
    [b, a] = butter(4,0.1/(fs/2),"high");
    vel_filtered = filtfilt(b,a,velocity);
    position = cumtrapz(t_zejl,vel_filtered);
    [b, a] = butter(4,0.5/(fs/2),"high");
    position= filtfilt(b,a,position);
    %pos_filtered = cumtrapz(t_zejl,pos_filtered);
    position = position-mean(position);
    [pos_peak, pos_peak_loc] = findpeaks(position);
    [~, pos_min_loc] = findpeaks(-position);
    pos_min = position(pos_min_loc);
    % figure()
    % plot(t_zejl,position,pos_peak_loc/fs,pos_peak,'O',pos_min_loc/fs,pos_min,'*')
    % if length(pos_min)~=length(pos_peak)
    %     if length(pos_min)>length(pos_peak)
    %         pos_min(end) = [];
    %     else
    %         pos_peak(end) = [];
    %     end
    % end

    %h_zejl(i) =   mean( abs(   max(pos_peak)-min(pos_min)  )); %mean( abs(pos_peak(2:2:end) - pos_peak(1:2:end))  ) ;
    index_peak_zejl = 1;

    while h_zejl(i)>0.10
        h_zejl(i) = mean( abs(   maxk(pos_peak,index_peak_zejl)-mink(pos_min,index_peak_zejl)  ));
        index_peak_zejl = index_peak_zejl +1;
    end

    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(namelist{i}(1:2), 'co')
        labcontrol.nsubject = idx_co;
        labcontrol.cadence(idx_co) = 60*fs/n_moe; %steps/minute
        labcontrol.meanstepLength(idx_co) = D_moe/M_moe; %meters/step
        labcontrol.walkingspeed(idx_co) = labcontrol.meanstepLength(idx_co)*fs/n_moe;
        mean_height = 1.64; %meters, given by the main paper
        leg_length = alpha*mean_height;
        step_length(i) = 2*sqrt(2*leg_length*h_zejl(i) - h_zejl(i)^2);
        labcontrol.gaitspeedZejlstra(idx_co) = step_length(i) * length(pos_peak)/t_zejl(end);



        sgl_withoutturns = sgl;
        sgl(idx_turns,:) = [];

        sgl(:,4:6) = [];
        
        [sgl(:,3),sgl(:,2),sgl(:,1)] = algo_Moe_Nilssen(sgl(:,3),sgl(:,2),sgl(:,1),'tiltAndNoG'); %TILT CORRECTION METHOD

        sgl(:,1) = (sgl(:,1) - mean(sgl(:,1)));
        sgl(:,2) = (sgl(:,2) - mean(sgl(:,2)));
        sgl(:,3) = (sgl(:,3) - mean(sgl(:,3)));


        sgl_medfilt = medfilt1(sgl,3);

        [b, a] = ellip(3, 0.01, 100, 0.25/(fs/2)); % fs is sampling frequency
        GA = filtfilt(b, a, sgl_medfilt); % zero-phase filtering to avoid phase distortion
        BA = sgl-GA;
        SMAc = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);
        labcontrol.SMA(idx_co) = mean(SMAc);



        % %%%%%%% PROCESSING AND BINARY SIGNAL OF FREQUENCY FILTER %%%%%%%%%
        % win_freqfilt = 5*fs;
        % binary_signal_freqfilt = NaN(1,m);
        % sgl_freqfilt = [zeros(win_freqfilt/2,1); sgl(:,1); zeros(win_freqfilt/2,1)];
        % L = length(sgl_freqfilt);
        % 
        % %%%%NO OVERLAP
        % energy = [];
        % counter = 1;
        % for h=1:win_freqfilt:L-win_freqfilt-1
        %     window = sgl_freqfilt(h:h+win_freqfilt-1);
        %     energy(counter) = bandpower(window,fs, [0.5 3])*(win_freqfilt/fs);
        %     counter = counter + 1;
        % end
        % meanenergy_c(idx_co) = mean(energy);



        idx_co = idx_co +1;

    else
        labfallers.nsubject = idx_fl;
        labfallers.cadence(idx_fl) = 60*fs/n_moe; %steps/minute
        labfallers.meanstepLength(idx_fl) = D_moe/M_moe; %meters/step
        labfallers.walkingspeed(idx_fl) = labfallers.meanstepLength(idx_fl)*fs/n_moe;
        mean_height = 1.61; %meters, given by the main paper
        leg_length = 0.5*mean_height;
        step_length(i) = 2*sqrt(2*leg_length*h_zejl(i) - h_zejl(i)^2);
        labfallers.gaitspeedZejlstra(idx_fl) = step_length(i) * length(pos_peak)/t_zejl(end);



        sgl_withoutturns = sgl;
        sgl(idx_turns,:) = [];

        sgl(:,4:6) = [];
        [sgl(:,3),sgl(:,2),sgl(:,1)] = algo_Moe_Nilssen(sgl(:,3),sgl(:,2),sgl(:,1),'tiltAndNoG'); %TILT CORRECTION METHOD

        sgl(:,1) = sgl(:,1) - mean(sgl(:,1));
        sgl(:,1) = sgl(:,1)/std(sgl(:,1));
        sgl(:,2) = sgl(:,2) - mean(sgl(:,2));
        sgl(:,2) = sgl(:,2)/std(sgl(:,2));
        sgl(:,3) = sgl(:,3) - mean(sgl(:,3));
        sgl(:,3) = sgl(:,3)/std(sgl(:,3));

        sgl_medfilt = medfilt1(sgl,3);

        [b, a] = ellip(3, 0.01, 100, 0.25/(fs/2)); % fs is sampling frequency
        GA = filtfilt(b, a, sgl_medfilt); % zero-phase filtering to avoid phase distortion
        BA = sgl-GA;
        SMAf = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);
        labfallers.SMA(idx_fl) = mean(SMAf);



        %         %%%%%%% PROCESSING AND BINARY SIGNAL OF FREQUENCY FILTER %%%%%%%%%
        % win_freqfilt = 5*fs;
        % binary_signal_freqfilt = NaN(1,m);
        % sgl_freqfilt = [zeros(win_freqfilt/2,1); sgl(:,1); zeros(win_freqfilt/2,1)];
        % L = length(sgl_freqfilt);
        % 
        % %%%%NO OVERLAP
        % energy = [];
        % counter = 1;
        % for h=1:win_freqfilt:L-win_freqfilt-1
        %     window = sgl_freqfilt(h:h+win_freqfilt-1);
        %     energy(counter) = bandpower(window,fs, [0.5 3])*(win_freqfilt/fs);
        %     counter = counter + 1;
        % end
        % meanenergy_f(idx_fl) = mean(energy);





        idx_fl = idx_fl +1;
    end


    display(i)


end


% Boxplots

% make a boxplot for x:
figure
boxplot(labcontrol.walkingspeed)
hold on
plot(ones(length(labcontrol.walkingspeed),1), labcontrol.walkingspeed,'O')
% make a boxplot for y, and specify that it belongs at 2 along the x-axis:
hold on
boxplot(labfallers.walkingspeed,'Positions',2)
hold on
plot(ones(length(labfallers.walkingspeed),1)*2, labfallers.walkingspeed,'O')
% update axes x-ticks and labels:
set(gca(),'XTick',[1 2],'XTickLabels',{'Control','Fallers'})


%% Outliers extraction

%control.walkingspeed(4) = [];

q1 = quantile(labcontrol.walkingspeed, 0.25);
q3 = quantile(labcontrol.walkingspeed, 0.75);
interquartile_range = q3 - q1;

llim = q1 - 1.5*interquartile_range; %lower limit
ulim = q3 + 1.5*interquartile_range; %upper limit

labcontrol.gaitspeed.outliers = (labcontrol.walkingspeed < llim | labcontrol.walkingspeed > ulim);
labcontrol.gaitspeed.outliers_idx = find(labcontrol.gaitspeed.outliers==1);

q1 = quantile(labfallers.walkingspeed, 0.25);
q3 = quantile(labfallers.walkingspeed, 0.75);
interquartile_range = q3-q1;

llim = q1 - 1.5*interquartile_range; %lower limit
ulim = q3 + 1.5*interquartile_range; %upper limit

labfallers.gaitspeed.outliers = (labfallers.walkingspeed < llim | labfallers.walkingspeed > ulim);
labfallers.gaitspeed.outliers_idx = find(labfallers.gaitspeed.outliers==1);


%% Result calculation

labfallers.gaitspeed.mean = mean(labfallers.walkingspeed);
labfallers.gaitspeed.std = std(labfallers.walkingspeed);

labcontrol.gaitspeed.mean = mean(labcontrol.walkingspeed);
labcontrol.gaitspeed.std = std(labcontrol.walkingspeed);

labfallers.gaitspeedZ.mean = mean(labfallers.gaitspeedZejlstra);
labfallers.gaitspeedZ.std = std(labfallers.gaitspeedZejlstra);

labcontrol.gaitspeedZ.mean = mean(labcontrol.gaitspeedZejlstra);
labcontrol.gaitspeedZ.std = std(labcontrol.gaitspeedZejlstra);

threshold_control = mean(labcontrol.SMA); %-0.5*std([control.SMA fallers.SMA]);
threshold_fallers = mean(labfallers.SMA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


labfallers.gaitspeed.mean_noout = mean(labfallers.walkingspeed(~labfallers.gaitspeed.outliers));
labfallers.gaitspeed.std_noout = std(labfallers.walkingspeed(~labfallers.gaitspeed.outliers));

labcontrol.gaitspeed.mean_noout = mean(labcontrol.walkingspeed(~labcontrol.gaitspeed.outliers));
labcontrol.gaitspeed.std_noout = std(labcontrol.walkingspeed(~labcontrol.gaitspeed.outliers));


labfallers.gaitspeedZ.mean_noout = mean(labfallers.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));
labfallers.gaitspeedZ.std_noout = std(labfallers.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));

labcontrol.gaitspeedZ.mean_noout = mean(labcontrol.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));
labcontrol.gaitspeedZ.std_out = std(labcontrol.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));
















% %% Subgoal
% 
% only_control = 1;
% 
% sgl = [];
% g = 9.81;
% fs = 100;
% win_sma =1*fs;  % 1-second window
% 
% folderpath = strcat(cd,'\long-term-movement-monitoring/');
% direc = dir(fullfile(folderpath,'*.dat*'));
% 
% idx_co = 1;
% idx_fl = 1;
% 
% 
% namelist = {};
% for i=1:length(direc)
%     namelist{i,1} = direc(i).name(1:length(direc(i).name)-4);
% end
% 
% 
% %threshold = threshold + 0.1;
% % 
% % namelist([51]) =  [];
% % namelist([end-3]) = [];
% 
% if only_control==1
%     start_analysis = 1;
%     end_analysis = 40;
% elseif only_control == 0
%         start_analysis = 41;
%         end_analysis = length(namelist);
% elseif only_control == 2
%         start_analysis = 1;
%         end_analysis = length(namelist);
% end
% 
% 
% for i=start_analysis:end_analysis
% tic
%     sgl = rdsamp(strcat('\long-term-movement-monitoring\',namelist{i}));
%     sgl = sgl *g;
%     [m,n] = size(sgl);
%     t = (0:m-1)/fs;
% 
% 
% 
%     % Filtering and preprocessing
%     % sgl= sgl-mean(sgl);
%     % sgl = detrend(sgl);
%     % [b, a] = butter(4,30/fs);
%     % sgl(:,1:3) = filtfilt(b,a,sgl(:,1:3));
% 
%     % figure()
%     % plot(t(1:24000), sgl(1:24000,3))
% 
%     %%%%%%% The following steps are described in 11-appendix-3 %%%%%%%%%%
% 
% 
%     sgl(:,4:6) = [];
%     [sgl(:,3),sgl(:,2),sgl(:,1)] = algo_Moe_Nilssen(sgl(:,3),sgl(:,2),sgl(:,1),'tiltAndNoG'); %TILT CORRECTION METHOD
% 
% 
%     sgl(:,1) = sgl(:,1) - mean(sgl(:,1));
%     sgl(:,1) = sgl(:,1)/std(sgl(:,1));
%     sgl(:,2) = sgl(:,2) - mean(sgl(:,2));
%     sgl(:,2) = sgl(:,2)/std(sgl(:,2));
%     sgl(:,3) = sgl(:,3) - mean(sgl(:,3));
%     sgl(:,3) = sgl(:,3)/std(sgl(:,3));
% 
% 
%     sgl_medfilt = medfilt1(sgl,3);
% 
%     [b, a] = ellip(3, 0.01, 100, 0.25/(fs/2)); % fs is sampling frequency
%     GA = filtfilt(b, a, sgl_medfilt); % zero-phase filtering to avoid phase distortion
% 
%     BA = sgl-GA;
% 
% 
%     %threshold = 0.153*g;
%     SMA = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);
%     %threshold = prctile(SMA,97.5);
%     %SMA = ( sqrt(sgl(:,1).^2 + sgl(:,2).^2 + sgl(:,3).^2) )/m;
% 
%     % figure
%     % subplot(2,1,1)
%     % plot(t,SMA)%t(1:6000),ones(1,6000)*control.threshold_SMA)
%     % ylabel('SMA')
%     % yline(threshold,'r','LineWidth',1)
%     % subplot(2,1,2)
%     % plot(t,BA(:,3))
%     % ylabel('BA AP component')
%     % yline(threshold,'r','LineWidth',1)
% 
% 
% 
%     %%%%%%% PROCESSING AND BINARY SIGNAL OF FREQUENCY FILTER %%%%%%%%%
%     win_freqfilt = 4*fs;
%     binary_signal_freqfilt = NaN(1,m);
%     %sgl_freqfilt = [zeros(win_freqfilt/2,1); sgl(:,3)/g; zeros(win_freqfilt/2,1)];
%     sgl_freqfilt = sgl(:,1); %%%%% or (:,3) ????
%     L = length(sgl_freqfilt);
% 
%     %%%%NO OVERLAP
%     c = 0;
%     h=1;
% 
%     while h+win_freqfilt <= L + 1
% 
%         window = sgl_freqfilt(h:h+win_freqfilt-1);
%         %[pxx, faxis] = pwelch(window, 2*fs ,[], 512, fs);
%         % figure()
%         % plot(faxis,pxx)
%         % xlim([0 6])
%         % energy_our = bandpower(pxx,faxis,[0.5 3],"psd")*(win_freqfilt/fs);
%         % % %
%         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %  freq_idx = (faxis >= 0.5) & (faxis <= 3);
%         %  energy_integral = trapz(faxis(freq_idx), pxx(freq_idx)) * win_freqfilt/fs;
%         %
%         % %Frequency resolution
%         % df = fs/1024;
%         %
%         % %Integrate the PSD over the band to get power (mean power)
%         % power_in_band = sum(pxx(freq_idx) * df);
%         %
%         % %Convert power to energy by multiplying by window duration (seconds)
%         %
%         % energy_in_band = power_in_band * win_freqfilt/fs;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         energy = bandpower(window, fs, [0.5 3])*(win_freqfilt/fs);
%         % 
%         % [b,a] = butter(4,[0.5/(fs/2) 3/(fs/2)],"bandpass");
%         % window_filtered = filtfilt(b,a,window);
%         % energy_new  = sum(window_filtered.^2);
% 
%         if energy > 0.05
%             binary_signal_freqfilt(h:h+win_freqfilt-1) = 1;
%         else
%             binary_signal_freqfilt(h:h+win_freqfilt-1) = 0;
%         end
%         %display(h)
%         h = h+floor(win_freqfilt); %%% /2 to obtain overlap
% 
%         if (h+win_freqfilt)>L && c==0
%             h = L-win_freqfilt+1;
%             c=1;
%         end
%     end
% 
%     %binary_signal_freqfilt = binary_signal_freqfilt(1:length(sgl)); %to improve
% 
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     %%%%%%%%%%%%%%%%%% BINARY SIGNAL SMA %%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     if strcmp(namelist{i}(1:2), 'CO')
%         threshold = threshold_control;
%     else
%         threshold = threshold_fallers;
%     end
% 
% 
%     actv_idx = find(SMA>threshold);
%     binary_signal = zeros(1,length(sgl));
%     binary_signal(actv_idx) = 1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     %%%%%% BINARY SIGNAL MERGING %%%%%%%%
%     binary_signal = double(logical(binary_signal) | logical(binary_signal_freqfilt));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
%     actv_duration_samples = {};        % Initialize result cell array
%     k = 1;                     % Start index for binary_signal
%     seg_counter = 1;           % Counter for activity segments
% 
% 
%     while k < length(binary_signal)
%         total_zero_counter = 0;
%         if binary_signal(k) == 1
%             count = 0;
%             count_zeros = 0;
%             j = k;
%             while (j < length(binary_signal)) && (binary_signal(j) == 1 || count_zeros<1*fs)
%                 count = count + 1;
%                 if binary_signal(j) == 1
%                     count_zeros = 0;
%                 end
%                 if binary_signal(j) == 0
%                     total_zero_counter = total_zero_counter + 1;
%                     count_zeros = count_zeros + 1;
%                 end
%                 j = j + 1;
%             end
%             if count > 60*fs %&& total_zero_counter<5*count/100
%                 idx_start = k;
%                 idx_end = j - 1;   % Since j points one past the last '1'
%                 actv_duration_samples{seg_counter} = [idx_start, idx_end];
%                 seg_counter = seg_counter + 1;
%             end
%             k = j;  % Skip to first index after segment
%         else
%             k = k + 1;    % Step forward if no activity
%         end
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%% AMBRA FEDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %     indices = [];
%     % for k=1:length(actv_duration_samples)
%     %     indices(k,1) =  actv_duration_samples{k}(1);
%     %     indices(k,2) = actv_duration_samples{k}(2);
%     % end
%     %
%     % steps=[];
%     % strides=[];
%     % weirdOutputs = [];
%     %
%     % for j = 1:size(indices,1)
%     %     acc_interval = BA(indices(j,1):indices(j,2),3);
%     %     [weirdOutput, N_steps, stride] = step_stride_duration(acc_interval,fs);
%     %     N_of_steps(j,1)  = N_steps;
%     %     strides(j,1)= stride;
%     %     weirdOutputs(j,1) = weirdOutput;
%     % end
% 
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
%     % toc
%     %
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % indexes = cell(1,1);
%     % reg = [];
%     % zeros = 0;
%     % totzero = 0;
%     % start = 1;
%     %
%     % for k=1:length(binary_signal)
%     %
%     %     if length(reg) >= 6000 && zeros == 100
%     %         indexes{end+1} = [start k-100];
%     %         reg = [];
%     %
%     %     end
%     %
%     %     if binary_signal(k) == 0
%     %         zeros = zeros +1;
%     %         totzero = totzero +1;
%     %     end
%     %
%     %     if binary_signal(k) == 1
%     %         if zeros == 100
%     %             start = k
%     %             reg = []
%     %         end
%     %         zeros = 0;
%     %     end
%     %     reg(end+1) = binary_signal(k);
%     %
%     % end
%     %
% 
% 
%     %%%%%%%%%%%%%%%%%
% 
%     [total_steps(i), steps_per_bout, cadence_spm, weirdoutput] = stepcount_MNandAvgStride(BA(:,3),fs,actv_duration_samples);
% 
%     actv_duration_time = NaN(1,length(actv_duration_samples));
%     for ii=1:length(actv_duration_samples)
%         actv_duration_time(ii) = diff(actv_duration_samples{ii})/fs;
%     end
% 
%     %average_stride_dur(i) = mean(actv_duration_time);
%     total_walking_time = sum(actv_duration_time); %in seconds
%     total_walking_percent = total_walking_time*100/t(end);
%     %display(total_walking_percent(i))
%     disp(i)
% 
%     % % % for kk=1:length(actv_duration_samples)-1
%     % % %     check(kk) = actv_duration_samples{kk+1}(1) - actv_duration_samples{kk}(2);
%     % % % end
%     % % % double_check{i} = find(check<0); %to check for superimposition
% 
% 
% 
% 
% 
%     if strcmp(namelist{i}(1:2), 'CO')
% 
% 
%         control.walking_percent(idx_co) = total_walking_percent;
%         control.total_steps(idx_co) = total_steps(i);
%         control.avgstride(idx_co) = weirdoutput;
%         %control_average_stride_duration(idx_co) = mean(strides);
%         %control_meanweirdoutp (idx_co) = mean(weirdOutputs);
%         %control_totnumsteps (idx_co) = sum(N_of_steps);
% 
% 
%         idx_co = idx_co +1;
% 
%     else
%         fallers.walking_percent(idx_fl) = total_walking_percent;
%         fallers.total_steps(idx_fl) = total_steps(i);
%         fallers.avgstride(idx_fl) = weirdoutput;
%         %fallers_average_stride_duration(idx_fl) = mean(strides);
%         %fallers_meanweirdoutp (idx_fl) = mean(weirdOutputs);
%         %fallers_totnumsteps (idx_fl) = sum(N_of_steps);
% 
% 
%         idx_fl = idx_fl +1;
%     end
% 
% toc
% end
% computing_time = toc;
% 
% if only_control == 1 || only_control == 2
%     control.meanwp = mean(control.walking_percent);
%     control.stdwp = std(control.walking_percent);
%     control.medianwp = median(control.walking_percent);
%     idx_outl_wp = findOutliers(control.walking_percent);
%     control.meanwpNOoutl = mean(control.walking_percent(idx_outl_wp));
%     control.stdwpNOoutl = std(control.walking_percent(idx_outl_wp));
% 
%     control.meants = mean(control.total_steps);
%     control.stdts = std(control.total_steps);
%     control.mediants = median(control.total_steps);
%     idx_outl_ts = findOutliers(control.total_steps);
%     control.meantsNOoutl = mean(control.total_steps(idx_outl_ts));
%     control.stdtsNOoutl = std(control.total_steps(idx_outl_ts));
% 
% 
%     control.meanavgstride = mean(control.avgstride);
%     control.stdavgstride = std(control.avgstride);
%     control.medianavgstride = median(control.avgstride);
%     idx_outl_avgstride = findOutliers(control.avgstride);
%     control.meanavgstrideNOoutl = mean(control.avgstride(idx_outl_avgstride));
%     control.stdavgstrideNOoutl = std(control.avgstride(idx_outl_avgstride));
% end
% 
% 
% if ~only_control || only_control == 2
%     fallers.meanwp = mean(fallers.walking_percent);
%     fallers.stdwp = std(fallers.walking_percent);
%     fallers.medianwp = median(fallers.walking_percent);
%     idx_outl_wp = findOutliers(fallers.walking_percent);
%     fallers.meanwpNOoutl = mean(fallers.walking_percent(idx_outl_wp));
%     fallers.stdwpNOoutl = std(fallers.walking_percent(idx_outl_wp));
% 
% 
%     fallers.meants = mean(fallers.total_steps);
%     fallers.stdts = std(fallers.total_steps);
%     fallers.mediants = median(fallers.total_steps);
%     idx_outl_ts = findOutliers(fallers.total_steps);
%     fallers.meantsNOoutl = mean(fallers.total_steps(idx_outl_ts));
%     fallers.stdtsNOoutl = std(fallers.total_steps(idx_outl_ts));
% 
% 
%     fallers.meanavgstride = mean(fallers.avgstride);
%     fallers.stdavgstride = std(fallers.avgstride);
%     fallers.medianavgstride = median(fallers.avgstride);
%     idx_outl_avgstride = findOutliers(fallers.avgstride);
%     fallers.meanavgstrideNOoutl = mean(fallers.avgstride(idx_outl_avgstride));
%     fallers.stdavgstrideNOoutl = std(fallers.avgstride(idx_outl_avgstride));
% end
% 
% 
% 
% 
% 
% 
% 
% if only_control == 1 || only_control == 2
% 
%     fprintf('\n===== RESULTS NON FALLERS =====\n');
%     fprintf('The mean value of the total walking percentage is %.2f %%\n', control.meanwpNOoutl);
%     fprintf('The standard deviation of the total walking percentage is %.2f \n', control.stdwpNOoutl);
%     fprintf('The mean value of the total steps is %.2f %%\n', control.meantsNOoutl);
%     fprintf('The standard deviation of the total steps is %.2f \n', control.stdtsNOoutl);
% 
% end
% 
% if ~only_control || only_control == 2
%     fprintf('\n===== RESULTS FALLERS =====\n');
%     fprintf('The mean value of the total walking percentage is %.2f %%\n', fallers.meanwpNOoutl);
%     fprintf('The standard deviation of the total walking percentage is %.2f \n', fallers.stdwpNOoutl);
%     fprintf('The mean value of the total steps is %.2f %%\n', fallers.meantsNOoutl);
%     fprintf('The standard deviation of the total steps is %.2f \n', fallers.stdtsNOoutl);
% end
% 
% fprintf('\n===== TOTAL COMPUTATIONAL TIME =====\n');
% fprintf('Total time elapsed: %.2f\n', computing_time);

