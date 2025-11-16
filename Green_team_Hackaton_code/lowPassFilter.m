function [acc] = lowPassFilter(acc,FS, cutoff_freq)
% Apply butterworth low pass filter to data
[C,D] = butter(1,cutoff_freq/(FS/2));
acc = filtfilt(C,D,acc); 
end

