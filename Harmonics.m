clear all
close all

%% Opens audio file
file = 'GBS_Project.wav';
[z,zfs]=audioread(file);
unsmooth_faxis =  zfs*(0:length(z)-1)/length(z);
num_peaks = 15; %Number of harmonics to remove from the signal, this needs to be specified by the user.
num_peaks_view = 9;
%% Find Peaks
ClearSignal=PeakRemover(file,num_peaks); %This is the time domain signal without peaks
ThreeP = AutoPeak(file); %This finds the first three peaks of the first three harmonics (including fundamental)
m=MidFinder(ThreeP,num_peaks);%This finds the midpoint between all the harmonics
m2=round(m(2)*zfs); %Finds the midpoint between the fundamental and the first harmonic in seconds
mend=round(m(end)*zfs); %Finds the last midpoint between harmonics in seconds

%% Filter and plots harmonics
ft_wins_hann = HarmonicFilt(z, zfs, num_peaks, m, unsmooth_faxis);

figure;
for i = 1:num_peaks_view
    semilogx(unsmooth_faxis, ft_wins_hann{i});
    hold on
end
label(1) = "Fundamental";
label(2:num_peaks_view) = char("H"+(2:num_peaks_view));
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
legend(label');
xlim([15 20000]);
grid on;