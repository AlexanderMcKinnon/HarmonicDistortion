function [THDdB, THD_perc, harm_perc, max_THD, freq_of_max_THD, min_THD, freq_of_min_THD] = calculateTHD(signal, signal_fs, windowedHarmonics)

final_harmonic = 10;
distortion_at = cell(40000, 10);
for frequency = 1:40000
    
    frequency_index = round(frequency.*length(signal)/signal_fs);   % find 
%   the index in the signal of each integer value of frequency from 1Hz to 
%   40kHz  

    distortion_at{frequency}(1) = 0;    % set the distortion due to the 
%   first harmonic to 0 for every frequency 

    for harmonic = 2:final_harmonic
       
        distortion_at{frequency}(harmonic) = 10^(windowedHarmonics{harmonic}(frequency_index)/20);
%       calculate the distortion due to each harmonic at every frequency 
%       Note: distortion is given in decibels
%       These are the arrays for the harmonic plots  
    end
end

THDdB = zeros(1, 40000);
for frequency = 1:40000
    d_squared{frequency} = distortion_at{frequency}.^2;
    num(frequency) = sqrt(sum(d_squared{frequency}));
    THDdB(frequency) = 20*log10(sum(distortion_at{frequency})); % Create an array of
%   the THD calculated at each frequency between 1Hz and 40kHz
end

% finding minimum THD
freq_of_min_THD = find(THDdB(15:20000) == min(THDdB(15:20000)));
min_THD = THDdB(freq_of_min_THD);

% finding maximum THD
freq_of_max_THD = find(THDdB(15:20000) == max(THDdB(15:20000)));
max_THD = THDdB(freq_of_max_THD);

%% THD calculated as percentage relative to dB level of fundamental 
for frequency = 1:40000
    frequency_index = round(frequency.*length(signal)/signal_fs);
    a = 10^(windowedHarmonics{1}(frequency_index)/20);
    b = 10^(THDdB(frequency)/20);
    THD_perc(frequency) = (10^(THDdB(frequency)/20))./(10^(windowedHarmonics{1}(frequency_index)/20));
%     THD_perc(frequency) = num(frequency)/(10^(windowedHarmonics{1}(frequency_index)/20));
end
for harmonic = 2:final_harmonic
    harm_perc{harmonic-1} = (10.^(windowedHarmonics{harmonic}./20))./(10.^(windowedHarmonics{1}./20));
end
end



