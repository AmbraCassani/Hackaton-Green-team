function actv_duration_time = samples_to_time(actv_duration_samples,fs)

actv_duration_time = NaN(1,length(actv_duration_samples));
for ii=1:length(actv_duration_samples)
    actv_duration_time(ii) = diff(actv_duration_samples{ii})/fs;
end

end