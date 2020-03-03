% code generated from the steps outlined on page 28 of
% https://pdfs.semanticscholar.org/bb2d/d8b44b8b90bbc5efbabd692e3bf71a7dd370.pdf 

[input, fs_file] = audioread('FemaleSpeech-16-8-mono-3secs.wav');
[impulseResponse, fs] = audioread('GBS_Project.wav');

% input signal passed through system 
output = conv(input, impulseResponse);

% time delay between the input and output is found, and the signals are
% time aligned 
tdelay = finddelay(input, output);
if tdelay > 0
    output(1:tdelay) = [];
elseif tdelay < 0
    tdelay = -1*tdelay;
    input(1:tdelay) = [];
end

% visually check that the signals are time aligned 
figure() 
plot(output)
hold on 
plot(input)


% input signal is analysed in frames of 30ms 
frames = fs_file*0.03;
sections = ceil(length(input)/frames);
splitInputSignal = cell(sections, 1);
counter = 0;
for i = 1:sections
    if (counter + frames) < length(input)
        splitInputSignal{i} = input(counter + 1:counter + frames);
        counter = frames;
    else 
        splitInputSignal{i}(:) = input(counter + 1:length(input));
    end
end

% output signal is analysed in frames of 30ms 
splitOutputSignal = cell(sections, 1);
counter = 0;
for i = 1:sections
    if (counter + frames) < length(input)
        splitOutputSignal{i} = output(counter + 1:counter + frames);
        counter = counter + frames;
    else 
        splitOutputSignal{i}(:) = output(counter + 1:length(input));
    end
end
n = 1323;
freqAxis_sound = fs*(0:frames-1)/frames;
freqAxis_IR = fs*(0:length(impulseResponse)-1)\length(impulseResponse);
test_freq = fs*(0:n - 1)/n;

dft_splitInput = cell(sections, 1);
dft_splitOutput = cell(sections, 1);
scaled_dft_splitOutput = cell(sections, 1);
freq_out = cell(sections, 1);
% A 1323 point DFT is performed over each frame i, and the relative peaks
% of the output signal are scaled to the input signal to remove offset/gain
% bias

figure()
% scale calculated after log and abs
for i = 1:sections
    dft_splitInput{i} = 20*log10(abs(fft(splitInputSignal{i}, n)));
    dft_splitOutput{i} = 20*log10(abs(fft(splitOutputSignal{i}, n)));
    max_input = max(dft_splitInput{i});
    max_output = max(dft_splitOutput{i});
    scale = max_input/max_output;
    scaled_dft_splitOutput{i} = dft_splitOutput{i}.*scale;
    [scaled_dft_splitOutput{i}, freq_out{i}] = rlogbarkranged(test_freq, scaled_dft_splitOutput{i}, 0, 95927);
end


