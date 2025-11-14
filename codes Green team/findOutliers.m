function [index_NOoutliers] = findOutliers(signal)
% The function find the indexes of the outliers of the signal

q1 = quantile(signal, 0.25);
q3 = quantile(signal, 0.75);
interquartile_range = q3 - q1;

llim = q1 - 1.5*interquartile_range; %lower limit
ulim = q3 + 1.5*interquartile_range; %upper limit

signal = (signal < llim | signal > ulim);
index_NOoutliers = find(signal==0);


end