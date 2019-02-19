function [time, f] = YIN_Estimator(signal, Fs, varargin)
%signal - input audio signal
%Fs - sampling rate
%checkSilent - if enabled - check for non-silent frames
%time,f - time vector and associated fundamental frequencies estimated

%windowsize - minimum frequency 40 Hz
window = round(0.025*Fs);
N = length(signal);
framesNum = ceil(N/window);
%zero pad signal to have enough frames
signal = [signal, zeros(1,window*framesNum - N)];
frames = zeros(framesNum, window);
start = 1;
%break into frames
for i = 1:framesNum
    frames(i,:) = signal(start:start + window - 1);
    start = start + window;
end

%step 1 - calculate difference function 
d = zeros(framesNum,window);
xTmp = [frames, zeros(framesNum,window)];
for tau = 0:window-1
    for j = 1:window  
         d(:,tau+1) = d(:,tau+1) + (xTmp(:,j) - xTmp(:,j+tau)).^2;         
    end
end

%step 2 - cumulative mean normalized difference function
d_norm = zeros(framesNum,window);
d_norm(:,1) = 1;

for i = 1:framesNum
    for tau = 1:window-1
        d_norm(i,tau+1) = d(i,tau+1)/((1/tau) * sum(d(i,1:tau+1)));
    end
end

%step 3 - absolute thresholding
lag = zeros(1,framesNum);
th = 0.1;
for i = 1:framesNum
    l = find(d_norm(i,:) < th,1);
    if(isempty(l) == 1)
        [v,l] = min(d_norm(i,:));
    end
    lag(i) = l;    
end

%step 4 - parabolic interpolation
period = zeros(1,framesNum);
time = zeros(framesNum,window);
f = zeros(framesNum,window);

for i = 1:framesNum
    if(lag(i) > 1 && lag(i) < window)
        alpha = d_norm(i,lag(i)-1);
        beta = d_norm(i,lag(i));
        gamma = d_norm(i,lag(i)+1);
        peak = 0.5*(alpha - gamma)/(alpha - 2*beta + gamma);
    else
        peak = 0;
    end
    period(i) = (lag(i)-1) + peak;
    f(i,:) = Fs/period(i)*ones(1,window);
    time(i,:) = ((i-1)*window:i*window-1)/Fs;    
end

if ~isempty(varargin)
    [f] = silent_frame_classification(frames, f);
end

f = reshape(f',1,framesNum*window);
time = reshape(time',1,framesNum*window);

end
