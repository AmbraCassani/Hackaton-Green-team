function [ degrees ] = getTurningDegreesDS( gyroV, FS )
% calculate turning degrees (in abs value)
%   GyroV is the angular velocity around the vertical axis ALREADY FILTERED in  degrees/second
winSamples=size(gyroV,1);
seconds=(1:winSamples)/FS;
degrees=abs(trapz(seconds,gyroV));
end

