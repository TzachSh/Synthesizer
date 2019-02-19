function [ f ] = silent_frame(frames, f)

framesNum = size(frames,1);
N = size(frames,2);

%this method works for noisy signals - find PSD of signal, and its spectral flatness
%the assumption is that silent frames have only noise and therefore a flat
%power spectrum

spectralFlatness = zeros(1,framesNum);
energy = zeros(1,framesNum);

for i = 1:framesNum
    [psdw, ~] = pwelch(frames(i,:)); 
    %this method works in case the signal isn't noisy.
    energy(i) = 20*log10(sum(frames(i,:).^2));
    %the spectral flatness method works only if signal is noisy
    spectralFlatness(i) = geomean(psdw)/mean(psdw);
end

normSpecFlatness = spectralFlatness./max(spectralFlatness);
threshold = 0.9;
silent = normSpecFlatness >= threshold | energy < -50;

s = find(silent == 1);
for i = 1:length(s)
    f(s(i),:) = zeros(1,N);
end

%make sure that partially silent frames are taken into account

for i = 2:framesNum-1
    %voiced frame preceded or followed by silent frame
    if(silent(i) == 0 && silent(i-1) == 1) 
        %returns signal indices with very low amplitude
        index = envelope_follower(frames(i,:),1);
        f(i,index) = zeros(1,length(index));
    elseif(silent(i) == 0 && silent(i+1) == 1)
        index = envelope_follower(frames(i,:),0);
        f(i,index) = zeros(1,length(index)); 
    end
end
    
end