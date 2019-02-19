function y = getADSR(A, D, S, R, Dur, Fs)  

    N = Dur * Fs;
    t = [0:N-1]/Fs;

    y = interp1([0 (0.01+A) (0.12+D) Dur - (0.601 - S) - R-0.02 Dur], [0 1 0.4 0.4 0], t); % !!      

end