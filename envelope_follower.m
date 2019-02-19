function [silentIndices] = envelope_follower(signal, indicate)

%amplitude envelope tracker for data x
%indicate - frame preceding or following silent frame
%1 - following(after), 0-preceding(before)
%Envelope Detection based on Hilbert Transform 

analys=hilbert(signal);
env=abs(analys);
N = length(signal);

%applying moving average filter to further smooth envelope

M = 10;
b = 1/M * ones(1,M);
a = 1;
env = filter(b,a,env);

%find signal indices with very low value of amplitude envelope
n = find(abs(env) < 0.5*max(abs(env)));
silentIndices = [];

if(~isempty(n))

    temp = [0 cumsum(diff(n)~=1)];
    elements = n(temp==mode(temp));

    %silent regions can only be at the beginning or end of frame, not in the
    %middle. The number of indices belonging to silent frames has to be
    %greater than 10

    if(numel(elements) >= 10)
        if(indicate == 1)
            silentIndices = 1:max(elements);
        else
            silentIndices = min(elements):N;
        end
    end

end

end