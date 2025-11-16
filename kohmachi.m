%an efficent way of smoothing the low frequnecy microtremor signals
%this code is compiled  from konno_ohmachi.py by Jian Shi 
%erman sent�rk&hliva
%  Original paper:
%         K. Konno & T. Ohmachi (1998) "Ground-motion characteristics estimated
%         from spectral ratio between horizontal and vertical components of
%         microtremor." Bulletin of the Seismological Society of America.
%         Vol.88, No.1, 228-241.
%[Inputs]
%         signal: Signal to be smoothed in frequency domain.
%         freq_array: Frequency array corresponding to the signal.
%                     It must have the same length as "raw_signal".
%         smooth_coeff: A parameter determining the degree of smoothing.
%                       The lower this parameter, the more the signal
%                       is smoothed.
%         
% 
%           (Note: "raw_signal" and "freq_array" can be Python lists, 1D numpy
%                  arrays, or 2D 1-column/1-row arrays. The data type affects the
%                  running time. For optimum speed, use 1D numpy array as input.)
%     [Output]
%         y: Smoothed signal (1D numpy array).
% The formula of Konno-Ohmachi smoothing window is here:
%         http://www.geopsy.org/wiki/index.php/Smoothing_details
% Copyright (c) 2013-2017, Jian Shi
function y = kohmachi(signal,freq_array,smooth_coeff)
    
    x = signal;
    f = freq_array;
    f_shifted = f/(1+1e-4);
    L = length(x);
    for i = 1 : L
        if (i ~= 1) && (i ~= L)
            
            z = f_shifted / f(i);
            w = ((sin(smooth_coeff * log10(z)) / smooth_coeff) ./ log10(z)) .^ 4;
            w(isnan(w)) = 0;
            y(i) = (w * x') / sum(w);
        end
        clearvars fc
    end
    y(1) = y(2);
    y(L) = y(L-1);
end