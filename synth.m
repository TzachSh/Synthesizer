function y=synth(app,freq,dur,amp,shift,Fs,type,VecADSR_Amp,ADSRFreq)
% y=synth(freq,dur,amp,Fs,type)
%
% Synthesize a single note
%
% Inputs:
%  freq - frequency in Hz
%  dur - duration in seconds
%  amp - Amplitude in range [0,1]
%  Fs -  sampling frequency in Hz
%  type - string to select synthesis type
%         current options: 'fm', 'sine', or 'saw'

% Copyright (c) 2009 Ken Schutte
% more info at: http://www.kenschutte.com/midi

if nargin<5
    error('Five arguments required for synth()');
end

N = floor(dur*Fs);

if N == 0
    warning('Note with zero duration.');
    y = [];
    return;
    
elseif N < 0
    warning('Note with negative duration. Skipping.');
    y = [];
    return;
end

n=0:N-1;
switch app.WaveType
    case 'Sine'
        if strcmp(app.FreqADSRSwitch.Value, 'On')
            envFreq = getADSR(app,app.VecADSR_Freq(1),app.VecADSR_Freq(2),app.VecADSR_Freq(3),app.VecADSR_Freq(4),dur,Fs);
            WaveOut = amp .* sin(2 * pi * freq/fs * envFreq .* n + shift); % fourier synthesis
        else
            WaveOut = amp .* sin(2 * pi * n * freq/Fs + shift); % fourier synthesis
        end
    case 'Sawtooth'
        if strcmp(app.FreqADSRSwitch.Value, 'On')
            envFreq = getADSR(app,app.VecADSR_Freq(1),app.VecADSR_Freq(2),app.VecADSR_Freq(3),app.VecADSR_Freq(4),dur,Fs);
            WaveOut = amp .* sawtooth(2 * pi * freq/Fs * envFreq .* n + shift);
        else
            WaveOut = amp .* sawtooth(2 * pi * n * freq/Fs + shift); % fourier synthesis
        end
    case 'Square'
        if strcmp(app.FreqADSRSwitch.Value, 'On')
            envFreq = getADSR(app,app.VecADSR_Freq(1),app.VecADSR_Freq(2),app.VecADSR_Freq(3),app.VecADSR_Freq(4),dur,Fs);
            WaveOut = amp .* square(2 * pi * freq/Fs * envFreq .* n + shift);
        else
            WaveOut = amp .* square(2 * pi * n * freq/Fs + shift); % fourier synthesis
        end
end

if strcmp(app.AmpADSRSwitch.Value, 'On')
    AmpEnv = getADSR(VecADSR_Amp(1),VecADSR_Amp(2),VecADSR_Amp(3),VecADSR_Amp(4),dur,Fs);
    y = AmpEnv .* WaveOut;
else
    y = WaveOut;
end
% smooth edges w/ 10ms ramp
if (dur > .02)
    L = 2*fix(.01*Fs)+1;  % L odd
    ramp = bartlett(L)';  % odd length
    L = ceil(L/2);
    y(1:L) = y(1:L) .* ramp(1:L);
    y(end-L+1:end) = y(end-L+1:end) .* ramp(end-L+1:end);
end