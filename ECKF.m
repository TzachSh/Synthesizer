function [f,amp,phase,xEst] = ECKF(signal,Fs,noiseParam, numFramesToWait)

%signal - incoming signal 
%Fs - sampling rate
%f - estimated pitch
%noiseParam - parameter for determining process noise
%numFramesToWait - number of buffers to wait after note onset - this will vary
%from instrument to instrument, depending on how strong its attack is
%amp - estimated amplitude of fundamental
%phi - estimated phase of fundamental
%xEst - estimated fundamental component of signal

%break signal into frames of 25 ms
frameLength = round(0.025*Fs);
framesNum = ceil(length(signal)/frameLength);
signal = [signal, zeros(1,framesNum*frameLength - length(signal))];

Ts = 1/Fs;
H =[0,0.5,0.5];

%Kalman filter variables
K = zeros(3,1);
flag = -1;
Kthres = 0.01;
xEst = zeros(1,length(signal));
amp = zeros(1,length(signal));
phase = zeros(1,length(signal));
Q = zeros(1,length(signal));
f = zeros(1,length(signal));
start = 1;
n = 1;
%state of current frame - initially silent
curFrame = 1;

while(start + frameLength - 1 < length(signal))
    yFrame = signal(start:start+frameLength-1);
    prevFrame = curFrame;
    curFrame = is_silent(yFrame);
    
    %if current frame is silent, then continue
    if(curFrame == 1)
        flag = 0;
        start = start + frameLength;
        n = start;
        continue;
        
    %transition from non-silent to silent frame
    elseif(prevFrame == 1 && curFrame == 0)
        %start counting buffers
        count = 0;
        %wait for a few frames before doing analysis to ignore 'attack'
        %of an instrument where frequencies go haywire
        while(count < numFramesToWait)
            count = count+1;
            start = start + frameLength;
        end
        if(start + frameLength - 1 < length(signal))
            yFrame = signal(start:start+frameLength-1);
            flag = 1;
            n = start;
        else
            break;
        end
    end
          
    %reset covariance matrix and calculate initial states
    if(flag == 1)

        %calculate initial state by taking an FFT and detecting the first peak
        n = n + 1;
        yBuf = yFrame;
        minf = length(yBuf);
        nfft = 2^nextpow2(4*minf);
        fBins = linspace(-Fs/2,Fs/2,nfft+1);
        win = blackman(minf);
        yBuf = (yBuf - mean(yBuf)).* win';
        Y = fftshift(fft(yBuf, nfft));
        %considering minimum possible frequency to be 50Hz, we ignore all
        %bins that are below 50Hz. Number of bins below 50Hz = 50/(fs/2*nfft)
        nbins_below50 = round(50/(Fs/2*nfft));
        Y = Y(nfft/2 + nbins_below50:end); 
        fBins = fBins(nfft/2 + 1 + nbins_below50:end);
        mag = abs(Y)./mean(win);
        phi = angle(Y);
        %take the least frequency peak to be fundamental
        [m, mPos] = findpeaks(mag);
        [val, index] = sort(m, 'descend');
        %we assume that fundamental frequency is the minimum of largest peaks frequency
        mPos = min(mPos(index(val > 0.5*max(val))));  
        [a1,ppos] = parabolic_interpolation(mag(mPos-1),mag(mPos),mag(mPos+1));
        a1 = a1/nfft;
        f1 = fBins(mPos) + (ppos * Fs/nfft);
        [phi1,pos] = parabolic_interpolation(phi(mPos-1),phi(mPos),phi(mPos+1));
        x0 = [exp(1i*2*pi*f1*Ts);a1*exp(1i*2*pi*f1*n*Ts + 1i*phi1);...
            a1*exp(-1i*2*pi*f1*n*Ts - 1i*phi1)];
        P0 = 0;        

        %reset covariance matrix
        if(abs(min(K)) < Kthres)
            P_last = P0;
            x_last = x0;
            flag = 0;
        end
    end
    
        for k = 1:frameLength
        %ekf equations
        K = (P_last*H')/(H*P_last*H' + 1);
        P = P_last - K*H*P_last;
        x = x_last + K*(yFrame(k) - H*x_last);
        x_next = [x(1);x(1)*x(2);x(3)/x(1)];
        F = [1,0,0;x(2),x(1),0;-x(3)/(x(1)^2),0,1/x(1)];
        
        %adaptive process noise based on error 
        Q(n) = 10^-(noiseParam-(abs(yFrame(k) - H*x)));
        
        P_next = F*P*F' + Q(n)*eye(3);
        f(n) = abs(log(x(1))/(1j*Ts*2*pi));
        amp(n) = abs(x(2));
        phase(n) = abs(-1i * (log(x(2)/amp(n))-(2*pi*f(n)*Ts*n)));
        xEst(n) = H*x;

        P_last = P_next;
        x_last = x_next;
        n = n + 1;        
        end
    start = start + frameLength;
end

end

