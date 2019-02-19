function[peak,pos] = parabolic_interpolation(alpha,beta,game)
    %given 3 points, it returns the result of their parabolic interpolation
    pos = 0.5 * ((alpha-game)/(alpha - 2*beta + game));
    peak = beta - 0.25*(alpha-game)*pos;
end