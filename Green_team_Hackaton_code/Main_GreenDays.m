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

namelist = create_list(fullfile(strcat(cd,'\long-term-movement-monitoring/LabWalks/'),'*.dat*'));


thresholdDegrees=45;%only turns with degrees equal or over this threshold are identified
thresholdLPfiltering=1.5;%Hz
thresholdPeakDS=15;%degrees/second
thresholdCrossingDS=5;%degrees/second
thresholdIntraTurnS=0.05;%seconds
thresholdTurnDurationLowS=0.5; %minimum turn duration
thresholdTurnDurationHighS=10; %maximum turn duration
fs = 100; %sampling frequency
win_sma = 1*fs; %window of the SMA computation
alpha = 0.5; % relation between leg length and height

[labcontrol, labfallers,threshold_control,threshold_fallers] = gait_speed(namelist,fs,win_sma,alpha,thresholdDegrees,thresholdLPfiltering,thresholdPeakDS,thresholdCrossingDS, thresholdIntraTurnS, thresholdTurnDurationLowS, thresholdTurnDurationHighS);




%% Subgoal

only_control = 2; % 1 = only control are analyzed, 0 = only fallers are analyzed, 2 : complete analysis (both control and fallers)

sgl = [];
g = 9.80665;
fs = 100; % sampling frequency
win_sma =1*fs;  % 1-second window
win_freqfilt = 4*fs; % window on which the energy computation is performed

control = struct();
fallers = struct();
idx_co = 1;
idx_fl = 1;

namelist = create_list(fullfile(strcat(cd,'\long-term-movement-monitoring/'),'*.dat*')); % create a cell containing a name of each file


[control, fallers] = subgoals(control, fallers, only_control, g, fs, win_sma, win_freqfilt, namelist, idx_co, idx_fl, threshold_control, threshold_fallers);

%% Classifier Logistic Regression with Lasso Regularization

log_regr('long-term-movement-monitoring\ClinicalDemogData_COFL.xlsx');