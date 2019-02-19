function [silent] = is_silent(frame)
%silent frame classification in music 

%method for clean signal - calculate signal energy and see if it is below
%a certain threshold

[psdw, w] = pwelch(frame); 
%this method works in case the signal isn't noisy.
energy = 20*log10(sum(frame.^2));
%the spectral flatness method works only if signal is noisy
spectralFlatness = geomean(psdw)/mean(psdw);

threshold = 0.9;
silent = spectralFlatness >= threshold | energy < -50;
end

