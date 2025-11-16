function binary_signal_freqfilt = compute_frequency_activity_mask(sgl, fs,win_freqfilt)

% Energy threshold filtering of the signal based on Vertical acceleration

    [m,~] = size(sgl);
    binary_signal_freqfilt = NaN(1,m);
  
    sgl_freqfilt = sgl(:,1); 
    L = length(sgl_freqfilt);

    %%%%NO OVERLAP
    c = 0;
    h=1;

    while h+win_freqfilt <= L + 1

        window = sgl_freqfilt(h:h+win_freqfilt-1);
       
        energy = bandpower(window, fs, [0.5 3])*(win_freqfilt/fs);
      

        if energy > 0.05
            binary_signal_freqfilt(h:h+win_freqfilt-1) = 1;
        else
            binary_signal_freqfilt(h:h+win_freqfilt-1) = 0;
        end
        
        h = h+floor(win_freqfilt); %%% /2 to obtain overlap

        if (h+win_freqfilt)>L && c==0 % this if is made to avoid of going of the signal length limits
            h = L-win_freqfilt+1;
            c=1;
        end
    end
end