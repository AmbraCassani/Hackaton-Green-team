function actv_duration_samples = extract_walking_bouts(binary_signal, fs)

 actv_duration_samples = {};        % Initialize result cell array
    k = 1;                     % Start index for binary_signal
    seg_counter = 1;           % Counter for activity segments


    while k < length(binary_signal)
        total_zero_counter = 0;
        if binary_signal(k) == 1
            count = 0;
            count_zeros = 0;
            j = k;
            while (j < length(binary_signal)) && (binary_signal(j) == 1 || count_zeros<1*fs)
                count = count + 1;
                if binary_signal(j) == 1
                    count_zeros = 0;
                end
                if binary_signal(j) == 0
                    total_zero_counter = total_zero_counter + 1;
                    count_zeros = count_zeros + 1;
                end
                j = j + 1;
            end
            if count > 60*fs %&& total_zero_counter<5*count/100 % this value can be used to avoid sections composed mainly by zeros
                idx_start = k;
                idx_end = j - 1;   % Since j points one past the last '1'
                actv_duration_samples{seg_counter} = [idx_start, idx_end];
                seg_counter = seg_counter + 1;
            end
            k = j;  % Skip to first index after segment
        else
            k = k + 1;    % Step forward if no activity
        end
    end


end