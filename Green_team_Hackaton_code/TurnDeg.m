function [ degr ] = TurnDeg( gyro_vert, fs )
% calculate turning degrees (in abs value)
%   GyroV is the angular velocity around the vertical axis ALREADY FILTERED in  degrees/second
L=size(gyro_vert,1);
t=(1:L)/fs;
degr=abs(   trapz(t,gyro_vert)   );
end

