clear all
close all
file = 'GBS_Project.wav';
[z,zfs]=audioread(file);
xaxis=transpose([0:1/zfs:(length(z)-1)/zfs]);
unsmooth_faxis =  zfs*(0:length(z)-1)/length(z);
num_peaks = 15;
last_harmonic = 6;  % Choose which harmonic is last to display - rest will be plotted but greyed out with noise


ClearSignal=PeakRemover(file,num_peaks);

ThreeP = AutoPeak(file);
m=MidFinder(ThreeP,num_peaks);

m2=round(m(2)*zfs);
mend=round(m(end)*zfs);

% figure(1);
% plot(xaxis,z,'r');
% hold on;
% plot(xaxis,ClearSignal,'b');
% xlabel('time(sec)');
% ylabel('AU magnitude');
% legend('Original singal', 'Signal with no harmonics');
% title('Signal impulse with and without harmonics removed');
%% Bark Log Smoothing

%[old_z, faxis] = rlogbark(unsmooth_faxis, old_z);
%z = interp1(faxis, old_z,unsmooth_faxis);

%% Create the Hanning windows for each harmonic
m_ind = round(m.*zfs);
for filt_num = [1:1:num_peaks]
   
    hann_temp = hann(m_ind(filt_num) - m_ind(filt_num + 1));
    c = zeros(m_ind(filt_num + 1), 1);
    d = zeros((length(z)-m_ind(filt_num)), 1);
    wins_hann{filt_num} = [c' hann_temp' d']'.*z;
end

%% Filter and FT hanning windows 

i = 1;
while i <= length(wins_hann)
    ft_wins_hann{i} = 20*log10(abs(fft(wins_hann{i})));
    [ft_wins_hann{i}, faxis] = rlogbark(unsmooth_faxis, ft_wins_hann{i});
    ft_wins_hann{i} = interp1(faxis, ft_wins_hann{i}, unsmooth_faxis);
    i = i + 1;
end




%% Plot the harmonics 
figure(3)
subplot(2, 1, 1)
% Plot the FT of the signal - with smoothing
semilogx(faxis, rlogbark(unsmooth_faxis, 20*log10(abs(fft(z)))),'LineWidth',2) ;
hold on
semilogx(unsmooth_faxis, 20*log10(abs(fft(z))));
xlim([15 20000]);
title('Plot of Fourier Transform of Impulse');
xlabel('Frequency (Hz)');
ylabel('Amplutide (dB)');

subplot(2, 1, 2)
for i = 1:9
    semilogx(unsmooth_faxis, ft_wins_hann{i});
    hold on
end
xlim([15 20000]);
% 
% 
% harmonics = sqrt(ft_wins_hann{2}.^2 + ft_wins_hann{3}.^2 + ft_wins_hann{4}.^2 + ft_wins_hann{5}.^2 + ft_wins_hann{6}.^2 + ft_wins_hann{7}.^2 + ft_wins_hann{8}.^2 + ft_wins_hann{9}.^2);
% figure(20)
% semilogx(faxis, harmonics);
% hold on 
% semilogx(faxis, sqrt(ft_wins_hann{1}.^2));
% title("Fourier Tranforms");
% xlabel("Frequency (Hz)");
% ylabel("Magnitude (dB)");
% legend("Square root of the sum of the squares of the harmonics (2 - 9)", "Magnitude of the fundamental harmonic");
% semilogx(faxis, (harmonics./sqrt(ft_wins_hann{1}.^2)));

%% THD
final_harmonic_thd = 10;
% 
% zFFT = fft(z);
% psd_z = periodogram(z, [], 'twosided');
% 
% % thds = cell(1, 3);
% thds_all = cell(1, 40000); 
% distortion_at = cell(1, 40000); 

%% THDs calculated using magnitude (dB) - with Christophe method


% for integer values of frequency between 1 and 40000Hz, get the magnitude
% of each harmonic at that frequency
% output is cell array: 1 column for each frequency value, each containing
% a vector of the magnitudes at that frequency  
for freqs = [1:40000]
    freq_ind = round(freqs.*length(z)/zfs);
    distortion_at{freqs}(1) = ft_wins_hann{1}(freq_ind);
    for harms = [2:final_harmonic_thd]
        temp = 10^((ft_wins_hann{harms}(freq_ind) - ft_wins_hann{1}(freq_ind))/20)*100;
        distortion_at{freqs} = [distortion_at{freqs} temp];
    end
end

for freqs = [1:40000]
    thd_mag_at(freqs) = sum(distortion_at{freqs});
end



% figure(3)
% subplot(2, 1, 1)
% % Plot the desired number of harmonics, and grey out the ones above to
% % represent noise 
% % j = 1;
% % while j <= length(ft_wins_hann)
% %     if j <= last_harmonic
% %         semilogx(faxis, ft_wins_hann{j});
% %         hold on
% %     else
% % %         semilogx(faxis, ft_wins_hann{j}, 'Color',[0,0,0] + 0.8, 'LineWidth', 0.5);
% %     end
% %     j = j + 1;
% % end
% % xlim([15 20000]);
% % title('Plot of Fourier Transform of Impulse Fundamental and Harmonics: Hanning Window');
% % xlabel('Frequency (Hz)');
% % ylabel('Amplutide (dB)');
% % subplot(2, 1, 2)
% % semilogx(thd_mag_at);
% % xlim([15 20000]);
% % title("THD by magnitudes")
% % xlabel("Frequency (Hz)");
% % ylabel("Distortion (%)");
% 
% % power of signal: abs(fft).^2/transform length
% powerFund = abs(ft_wins_hann{1}).^2/length(ft_wins_hann{1});
% energyFund = powerFund.*unsmooth_faxis;
% 
% powerHarms = cell(1, final_harmonic_thd - 1);
% energyHarms = cell(1, final_harmonic_thd - 1);
% 
% for j = [2:final_harmonic_thd]
%     powerHarms{j - 1} = abs(ft_wins_hann{j - 1}).^2/length(ft_wins_hann{j - 1});
%     energyHarms{j - 1} = powerHarms{j - 1}.*unsmooth_faxis;
%     energyHarms{j - 1} = energyHarms{j - 1}.^2;
% end
% 
% harmonicEnergy = [];
% 
% 
% for freq = [1:1:length(faxis)]
%     sum = 0;
%     for j = [1:1:final_harmonic_thd - 1]
% %         temp = energyHarms{j}(freq);
%         sum = energyHarms{j}(freq) + sum;
%     end
%     harmonicEnergy(freq) = sum;
% end
% 
% numerators = sqrt(harmonicEnergy);
% thd_byEnergy = [];
% 
% thd_byEnergy = (numerators./energyFund').*100;
% 
% figure(1)
% semilogx(faxis, thd_byEnergy);
% % 
% % 
% % for freqs = [1:40000]
% %     freq_ind = round(freqs.*length(z)/zfs);
% %     mag_of_harmonic_at{freqs}(1) = ft_wins_hann{1}(freq_ind);
% %     for harms = [2:final_harmonic_thd]
% %         temp = (ft_wins_hann{harms}(freq_ind));
% %         mag_of_harmonic_at{freqs} = [mag_of_harmonic_at{freqs} temp];
% %     end
% % end
% 
% % Create the numerator of the calculation of the THD 
% % Creates a numerator for each frequency 
% % numers_test_mag = [];
% % for f = [1:40000]
% %     sum_test = 0;
% %     for l = (2:1:final_harmonic_thd)
% %         sum_test = sum_test + (mag_of_harmonic_at{f}(l))^2;
% %     end
% %     numers_test_mag(f) = sum_test;
% % end
% % 
% % thd_all_mag = zeros(1, 40000);
% % % Calculate THD-F
% % for f = [1:40000]
% %     thd_all_mag(f) = sqrt(numers_test_mag(f))/mag_of_harmonic_at{f}(1);
% % end
% % 
% % thd_all_mag = thd_all_mag*100;
% 
% figure(2)
% semilogx(thd_mag_at);
% hold on 
% % semilogx(thd_all_mag);
% xlim([15 20000]);

