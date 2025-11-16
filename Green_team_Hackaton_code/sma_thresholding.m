function binary_signal_SMA = sma_thresholding(sgl,threshold_control,threshold_fallers,SMA,namelist,i)

% The following function takes the threshold value previously computed for
% the gait speed analysis and use it to find activity period

if strcmp(namelist{i}(1:2), 'CO')
    threshold = threshold_control;
else
    threshold = threshold_fallers;
end


actv_idx = find(SMA>threshold);
binary_signal_SMA = zeros(1,length(sgl));
binary_signal_SMA(actv_idx) = 1;

end