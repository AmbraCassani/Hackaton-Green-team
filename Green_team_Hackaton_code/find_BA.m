function [BA, sgl] = find_BA(sgl,idx_turns,fs,NoTurns)

% Discern Body activity component from gravitational part

if NoTurns && ~isempty(idx_turns)
    sgl(idx_turns,:) = [];
end

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


end