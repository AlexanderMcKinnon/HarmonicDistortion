%% Non-App Version
clear all 
close all 
% Currently hardcoded for input file GBSTest.wav
% function Harmonics(callmode)
    %% Opens audio file
    
    file = 'GBS_Project.wav';   % Original IR

%     setpath;
%     current = AppConversion;
%     cd(MEASURES_SWEEP_PATH);
% 
%     if strcmpi(callmode,'main')
%         [file,path] = uigetfile('*.wav','Select a wav file to load');
% %         audio = strcat(path,file);
% 
%     elseif strcmpi(callmode,'plotfig')
%         file = evalin('base','IR_file');
%         path = evalin('base','IR_path');
% %         audio = strcat(path,file);
%     end
%     cd(current)
%     file = 'GBSTest.wav';
%     if not (file == 0) 
        [z,zfs]=audioread(file);
        % z = [zeros(1e5, 1); z; zeros(1e7, 1)];    % padding added to check if
        % works with different length signals
        unsmooth_faxis =  zfs*(0:length(z)-1)/length(z);
        backgroundNoiseFile = z(110000:end); %In the app this will be a seperate file that will have to be opened with audrioread. Christophe wants to record a signal to use as background noise.
        %% Variables to be in the User Interface
        num_peaks = 15; %Number of harmonics to remove from the signal, this needs to be specified by the user.
        num_peaks_view = 10;%Number of harmonics to view individually
        num_peaks_view_noise = 0; %number of harmonics to view as noise
        highlight_harmonic = 1; %Variable that determines which harmonic to make bold, if 0 then higlight none

        if (num_peaks_view > num_peaks); num_peaks_view = num_peaks; end
        if (num_peaks_view_noise > num_peaks); num_peaks_view_noise = num_peaks; end
        if (highlight_harmonic > num_peaks); highlight_harmonic = 0; end

        %% Remove noise
        [noise, noiseless_z] = removeNoise (backgroundNoiseFile, z, zfs, unsmooth_faxis);

        %% Find Peaks
        ClearSignal=PeakRemover(noiseless_z,zfs,num_peaks); %This is the time domain signal without peaks
        ThreeP = AutoPeak(noiseless_z,zfs); %This finds the first three peaks of the first three harmonics (including fundamental)
        m=MidFinder(ThreeP,num_peaks);%This finds the midpoint between all the harmonics
        m2=round(m(2)*zfs); %Finds the midpoint between the fundamental and the first harmonic in seconds
        mend=round(m(end)*zfs); %Finds the last midpoint between harmonics in seconds

        %% Filter and plots harmonics

        ft_wins_hann = HarmonicFilt(noiseless_z, zfs, num_peaks, m, unsmooth_faxis);% Function returns all the harmonics that have been windowed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ft_wins_hann_noise = 0; %Plots noise harmonics
        if num_peaks_view_noise > num_peaks_view
            for i = num_peaks_view+1:num_peaks_view_noise
                ft_wins_hann_noise = ft_wins_hann_noise + 10.^(ft_wins_hann{i}./20);
            end
        end
        ft_wins_hann_noise = ft_wins_hann_noise + 10.^(noise./20);
        ft_wins_hann_noise = 20*log10(ft_wins_hann_noise);
        
        %% THD
        [THDdB, THDperc, harm_perc, max_THD, freq_of_max_THD, min_THD, freq_of_min_THD] = calculateTHD(noiseless_z, zfs, ft_wins_hann);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         assignin('base','THD',THD);
%         assignin('base','ft_wins_hann_noise',ft_wins_hann_noise);
%         assignin('base','num_peaks_view',num_peaks_view);
%         assignin('base','highlight_harmonic',highlight_harmonic);
%         assignin('base','unsmooth_faxis',unsmooth_faxis);
%         assignin('base','ft_wins_hann',ft_wins_hann);
% 
%         if strcmpi(callmode,'main')
%             hold on
            figure;%Plots all the harmonics that is determined by user (num_peaks_view)
%             for i = 1:num_peaks_view
%                 if i == highlight_harmonic 
%                     semilogx(unsmooth_faxis, ft_wins_hann{i},'LineWidth',1.5); %This plots the Emphasised harmonic
%                 else
%                    semilogx(unsmooth_faxis, ft_wins_hann{i}); %This plots the rest of the harmonics
%                 end
%                 hold on
%             end
%             label(1) = "Fundamental"; %Creates a legend for graph
%             label(2:num_peaks_view) = char("H"+(2:num_peaks_view));
%             
            semilogx(THDdB);
%             label(length(label)+1)="THD (dB)";
            
%             semilogx(unsmooth_faxis, ft_wins_hann_noise,'LineStyle','--'); %'Color','#808080',
%             label(num_peaks_view+1) = "THD";
            
            title(['Harmonics for file ', file]); %Edits graph
            xlabel('Frequency (Hz)');
            ylabel('Amplitude (dB)');
%             legend(label');
            legend("THD (dB)")
            xlim([15 40000]);
            grid on;
            hold off
%         end



figure;%Plots all the harmonics that is determined by user (num_peaks_view)
% semilogx([1:1:40000],  100*THDperc, 'LineWidth', 1.5);
h = axes;
hold on
count = 0;
for i = 2:2:num_peaks_view
   count = count + 1;
   semilogx(unsmooth_faxis, 100*harm_perc{i-1}); %This plots the rest of the harmonics
   hold on 
end
%Create a legend for graph
% label2(1) = "THD (%)";

%             semilogx(unsmooth_faxis, ft_wins_hann_noise,'LineStyle','--'); %'Color','#808080',
%             label(num_peaks_view+1) = "THD";
% label2(1:count) = label(3:2:num_peaks_view);
% label2 = txt;
hold on 
semilogx(unsmooth_faxis, ones(1,length(unsmooth_faxis)),'LineStyle','--', 'Color', 'r');
title(['Harmonics % for even harmonics: file ', file]); %Edits graph
xlabel('Frequency (Hz)');
ylabel('Percentage with respect to H1 (%)');
% legend(label2');
legend('H2', 'H4', 'H6', 'H8', 'H10', '1% Threshold');
xlim([15 20000]);
set(h,'xscale','log')
grid on;
hold off
figure;%Plots all the harmonics that is determined by user (num_peaks_view)
% semilogx([1:1:40000],  100*THDperc, 'LineWidth', 1.5);
hold on
count = 0;
for i = 3:2:num_peaks_view
   count = count + 1;
   semilogx(unsmooth_faxis, 100*harm_perc{i-1}); %This plots the rest of the harmonics
   hold on
%    label3(count) = char("H"+(i));
end
%Create a legend for graph
% label2(1) = "THD (%)";
% label3(1:count) = char("H"+(3:2:num_peaks_view));

%             semilogx(unsmooth_faxis, ft_wins_hann_noise,'LineStyle','--'); %'Color','#808080',
%             label(num_peaks_view+1) = "THD";
hold on 
semilogx(unsmooth_faxis, 0.1*ones(1,length(unsmooth_faxis)),'LineStyle','--', 'Color', 'r');
title(['Harmonics % for odd harmonics: file ', file]); %Edits graph
xlabel('Frequency (Hz)');
ylabel('Percentage with respect to H1 (%)');
% legend(label3');
legend('H3', 'H5', 'H7', 'H9', '0.1% Threshold');
xlim([15 20000]);
set(h,'xscale','log')
grid on;
hold off
% end

%% set everything below the thresholds to 0% - assuming imperceptible 
for harmonic = [2:num_peaks_view]
    if rem(harmonic, 2) == 0    % even harmonic
        harm_perc{harmonic-1}(harm_perc{harmonic-1}<0.01)=0.00000001;
    else
        harm_perc{harmonic-1}(harm_perc{harmonic-1}<0.001)=0.00000001;
    end
end

figure;%Plots all the harmonics that is determined by user (num_peaks_view)
% semilogx([1:1:40000],  100*THDperc, 'LineWidth', 1.5);
h = axes;
hold on
count = 0;
for i = 2:2:num_peaks_view
   count = count + 1;
   semilogx(unsmooth_faxis, 100*harm_perc{i-1}); %This plots the rest of the harmonics
   hold on 
end
%Create a legend for graph
% label2(1) = "THD (%)";

%             semilogx(unsmooth_faxis, ft_wins_hann_noise,'LineStyle','--'); %'Color','#808080',
%             label(num_peaks_view+1) = "THD";
% label2(1:count) = label(3:2:num_peaks_view);
% label2 = txt;
hold on 
semilogx(unsmooth_faxis, ones(1,length(unsmooth_faxis)),'LineStyle','--', 'Color', 'r');
title(['Harmonics % for even harmonics: file ', file]); %Edits graph
xlabel('Frequency (Hz)');
ylabel('Percentage with respect to H1 (%)');
% legend(label2');
legend('H2', 'H4', 'H6', 'H8', 'H10', '1% Threshold');
xlim([15 20000]);
set(h,'xscale','log')
grid on;
hold off
figure;%Plots all the harmonics that is determined by user (num_peaks_view)
% semilogx([1:1:40000],  100*THDperc, 'LineWidth', 1.5);
h = axes;
hold on
count = 0;
for i = 3:2:num_peaks_view
   count = count + 1;
   if i == 3
        semilogx(unsmooth_faxis, 100*harm_perc{i-1}, 'LineWidth', 1.5);
   end
   semilogx(unsmooth_faxis, 100*harm_perc{i-1}); %This plots the rest of the harmonics
   hold on
%    label3(count) = char("H"+(i));
end
%Create a legend for graph
% label2(1) = "THD (%)";
% label3(1:count) = char("H"+(3:2:num_peaks_view));

%             semilogx(unsmooth_faxis, ft_wins_hann_noise,'LineStyle','--'); %'Color','#808080',
%             label(num_peaks_view+1) = "THD";
hold on 
semilogx(unsmooth_faxis, 0.1*ones(1,length(unsmooth_faxis)),'LineStyle','--', 'Color', 'r');
title(['Harmonics % for odd harmonics: file ', file]); %Edits graph
xlabel('Frequency (Hz)');
ylabel('Percentage with respect to H1 (%)');
% legend(label3');
legend('H3', 'H5', 'H7', 'H9', '0.1% Threshold');
xlim([15 20000]);
set(h,'xscale','log')
grid on;
hold off
% end

for i = 1:(num_peaks_view-1)
    a = 20*log10(harm_perc{i});
    percieved_harm_dB{i} = ft_wins_hann{1} + a;
end

figure;%Plots all the harmonics that is determined by user (num_peaks_view)
for i = 1:num_peaks_view
    if i == 1
        semilogx(unsmooth_faxis, ft_wins_hann{i}); %This plots the rest of the harmonics
        hold on
    elseif i == 3
        semilogx(unsmooth_faxis, percieved_harm_dB{i-1}, 'LineWidth', 1.5);
    else
        semilogx(unsmooth_faxis, percieved_harm_dB{i-1}); %This plots the rest of the harmonics
        hold on
    end
end
labelx(1) = "Fundamental"; %Creates a legend for graph
labelx(2:num_peaks_view) = char("H"+(2:num_peaks_view));

title(['Harmonics: adjusted for % thresholds']); %Edits graph
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
legend(labelx');
xlim([15 40000]);
grid on;
hold off

%% THD in dB from adjusted harmonic thresholds 
for i = 1:num_peaks_view
    if i == 1
        adjusted_harmonics{i} = ft_wins_hann{i};
    else 
        adjusted_harmonics{i} = percieved_harm_dB{i-1};
    end
end
[newTHDdB, newTHDperc, ~, ~, ~, ~, ~] = calculateTHD(noiseless_z, zfs, adjusted_harmonics);
figure;%Plots all the harmonics that is determined by user (num_peaks_view)
semilogx(newTHDdB);
title(['THD (dB): harmonics adjusted for % thresholds']); %Edits graph
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
legend("THD (dB)");
xlim([15 20000]);
grid on;
hold off


figure;%Plots all the harmonics that is determined by user (num_peaks_view)
semilogx(100*newTHDperc);
title(['THD (dB): harmonics adjusted for % thresholds']); %Edits graph
xlabel('Frequency (Hz)');
ylabel('%');
legend("THD (%)");
xlim([15 20000]);
grid on;
hold off

%     else
%         fprintf('Operation cancelled.\n');
%     end
%     cd(current)
