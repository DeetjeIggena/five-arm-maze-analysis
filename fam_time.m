% calculates time/ trial duration by substracting time-stamps, 
% alternative: use sampling rate & data length

function [duration, time ]= fam_time(a,b, data_length, samplingRate)

% calculates the time from the beginning to the end of the trial
duration = b-a;

time = samplingRate * data_length;

end





