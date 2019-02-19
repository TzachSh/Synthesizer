function detectedFreq = ZeroCrossing(signalOutput, Fs)
    % Count the number of sign changes between each adjencies.
    % The calculation requires only a half cycle.
    ZCR = sum(abs(diff(sign(signalOutput(:)))))/(2*length(signalOutput));
    F0_ZC=(ZCR * length(signalOutput) * Fs) / (2*length(signalOutput));
    detectedFreq = round(F0_ZC);
end