function [labcontrol,labfallers,threshold_control,threshold_fallers] = gait_speed(namelist,fs,win_sma,alpha,thresholdDegrees,thresholdLPfiltering,thresholdPeakDS,thresholdCrossingDS, thresholdIntraTurnS, thresholdTurnDurationLowS, thresholdTurnDurationHighS)

% The function has the main goal of computing the gait speed using Zejlstra
% method

h_zejl=ones(length(namelist),1);
labcontrol = struct();
labfallers = struct();
idx_co = 1;
idx_fl = 1;


for i= 1:length(namelist)
    %sgl = rdsamp('ltmm/CO001'); %use this to download it from the remote server

    sgl = rdsamp(strcat('long-term-movement-monitoring/LabWalks/',namelist{i}));
    g = 9.81; %needed for tilt correction method
    sgl(:,1:3) = sgl(:,1:3) * g; %convert to m/s^2
    [m,n] = size(sgl);

    t = (0:m-1)/fs;

    %% Turning detection to improve SMA computation

    [idx_turns] = findTurns(sgl(:,4), fs, thresholdDegrees, ...
        thresholdLPfiltering, thresholdPeakDS,thresholdCrossingDS,thresholdIntraTurnS, ...
        thresholdTurnDurationLowS, thresholdTurnDurationHighS );
    %%%% END OF TURNING DETECTION


    %%%%%%%%%%% Zejlstra gait speed method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [h_zejl,t_zejl,pos_peak] = find_hzejil(h_zejl,sgl,i,idx_turns,fs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    [labcontrol,labfallers,idx_co,idx_fl] = memorization(labcontrol,labfallers,namelist,alpha,h_zejl,...
        t_zejl,i,idx_co,idx_fl,sgl,idx_turns,fs,pos_peak,win_sma);


    display(i)


end


plot_Boxplot(labcontrol,labfallers,'gaitspeedZejlstra')



%% Result calculation


labfallers.gaitspeedZ.mean = mean(labfallers.gaitspeedZejlstra);
labfallers.gaitspeedZ.std = std(labfallers.gaitspeedZejlstra);

labcontrol.gaitspeedZ.mean = mean(labcontrol.gaitspeedZejlstra);
labcontrol.gaitspeedZ.std = std(labcontrol.gaitspeedZejlstra);
%%%%%%%%%%%%%%% MEAN SMA computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold_control = mean(labcontrol.SMA);
threshold_fallers = mean(labfallers.SMA);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%% Outliers removal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labfallers.gaitspeedZ.mean_noout = mean(labfallers.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));
labfallers.gaitspeedZ.std_noout = std(labfallers.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));

labcontrol.gaitspeedZ.mean_noout = mean(labcontrol.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));
labcontrol.gaitspeedZ.std_out = std(labcontrol.gaitspeedZejlstra(findOutliers(labfallers.gaitspeedZejlstra)));

end