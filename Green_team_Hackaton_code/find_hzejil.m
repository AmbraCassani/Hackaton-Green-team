function [h_zejl,t_zejl,pos_peak] = find_hzejil(h_zejl,sgl,i,idx_turns,fs)

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
[pos_peak, ~] = findpeaks(position);
[~, pos_min_loc] = findpeaks(-position);
pos_min = position(pos_min_loc);
index_peak_zejl = 1;

while h_zejl(i)>0.10
    h_zejl(i) = mean( abs(   maxk(pos_peak,index_peak_zejl)-mink(pos_min,index_peak_zejl)  ));
    index_peak_zejl = index_peak_zejl +1;
end


end