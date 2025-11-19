
clc; clear; close all; %clc: clears window, clear: clears all variabls, cloase all: closes alll figures
fc = 10000;      %Carrier frequency in Hertz (10khz)
fs = 100000;     %Sampling frequency in Hertz (100khz) , must be at least double the highest fequency 
M = 4;        % QPSK = 4 , " 4 possiable phases: 0, 90, 180, 270" 
bitspersymbol = log2(M); %defines how many bits of information are encoded into each transmitted symbol, 4 symbols, each carries 2 bits
symbolrate = 1000;   % defines how many symbols are transmetted per second
Ts = 1/symbolrate;   % duration of each symbol
% Binary-to-Decimal Function
BinarytoDecimal = @(binMat) sum(binMat .* 2.^(bitspersymbol-1:-1:0), 2); % converts binery to decimal , 2.^(bitspersymbol-1:-1:0) 2^1, 2^0 for weighting bits 
%Example: [1 0] becomes 1×2 + 0×1 = 2 , ( used this function instead of
%built in matlab function since matlab online does not support
%commuincation toolbox )

% Creating a randomw signal
numofbits = 100; %generates a 100 random bits
data = randi([0 1], numofbits, 1);  % Generates a 100×1 column vector of random integers (0 or 1)
chunkedBits = reshape(data, bitspersymbol, [])'; %organizes 100 bits into a 2×50 matrix (2 bits per column)
datasymbol = BinarytoDecimal(chunkedBits); % Converts each 2-bit row to a decimal symbol (0, 1, 2, or 3)

% PSK Modulation
% --------------------------------------------------------------
modData = exp(1j * 2*pi * (datasymbol) / M); %exp(1j * theta) is euler's formula = cos(theta) + j·sin(theta) "creats complex numbers"
%2pi × datasymbol / M: maps symbols {0,1,2,3} to phases {0, 90, 180, 270}
%results in a complex constellation points on the unit circle

% Time vector for modulated signal
t = 0:1/fs:length(modData)*Ts-1/fs; % generats a time vector from 0 to total transmission time, sampled at fs
transSignal = real(repelem(modData, fs*Ts) .* exp(1j*2*pi*fc*t).'); %repeats each symbol 100 times (fs×Ts = 100000×0.001) to match sampling rate
% exp(1*j*2*pi*fc*t): Generates carrier wave at 10kHz
% real() extracts the real part and leaves the complex part​

%manual implenmintation of awgn (matlab online does not support
%communication tool box)
SNR = 20; % signal to noise ratio
signal_power = mean(abs(transSignal).^2); %calculates signal power
noise_power = signal_power / (10^(SNR/10)); %calculates noise power
noise = sqrt(noise_power) * randn(size(transSignal)); %convirts snr from db to linear form
received_signal = transSignal + noise;  %received signal = transmitted signal + noise

% PSK Demodulation
numberofsymbols = length(modData); %total symbols to demodulate 
rxsymbols = zeros(numberofsymbols,1); % creates an array of size eqaul to "numberofsymbols"
refPhases = exp(1j * 2 * pi * (0:M-1)/M); % reference constellation points for all 4 QPSK phases

for i = 1:numberofsymbols    %loops through all of the recived symbols
    segment = received_signal((i-1)*fs*Ts+1:i*fs*Ts); %
    [~, val] = min(abs(mean(segment) - refPhases));   % extract symbols from received signal
    rxsymbols(i) = val - 1; % for 0-based indexing    %
end
%symbol by symbol demodulation 

rxBitsMat = rem(floor(rxsymbols(:) * ones(1, bitspersymbol) ./ 2.^(bitspersymbol-1:-1:0)),2);   %
rxBits = rxBitsMat';                                                                            % converts symbols to binary
rxBits = rxBits(:);                                                                             %

% plotting the random signal
figure;
subplot(3,1,1);
stem(data(1:20),'LineWidth',1.5); title('Original Binary Data (1st 20 bits)');
subplot(3,1,2);
plot(t, transSignal); title('PSK Modulated Random Signal');
subplot(3,1,3);
stem(rxBits(1:20),'LineWidth',1.5); title('Demodulated Binary Data (1st 20 bits)');

% plot in frequency domain 
txFreq = fftshift(fft(transSignal));
f = (-length(transSignal)/2:length(transSignal)/2-1)*(fs/length(transSignal));
figure;
plot(f, abs(txFreq));
title('Frequency Spectrum of PSK Modulated Signal');

%PSK Constellation Plot
figure;
plot(real(modData), imag(modData), '*');
axis equal;
grid on;
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('PSK Constellation Diagram (Random Data)');
xlabel('In-Phase'); ylabel('Quadrature');











% Voice signal
[voice, fsVoice] = audioread('RecordingPOC.wav');   % audioread: reads audion from file in between ' ', 
voice = resample(voice, fs, fsVoice);      %resample : converts the smapling frequency from the voice message to 100Khz (fs) , sampling frequency
disp('Playing The Original Voice Message...');  
sound(voice, fs); 
pause(length(voice)/fs + 1);


voiceNorm = voice / max(abs(voice)); %normilize amplitude to -1 to 1
voiceBinary = round((voiceNorm + 1) * 127.5);      %quantization to 8 bits (0-255)
voiceBitsMat = rem(floor(voiceBinary(:) * ones(1,8) ./ 2.^(7:-1:0)),2); % converts 8 bits decimal to binary
voiceBits = voiceBitsMat'; % converts 8 bits decimal to binary
voiceBits = voiceBits(:);% converts 8 bits decimal to binary

% Zero padding (adding extra zeros to a data array or at the end to meet a
% achive a desired length)
totalBits = length(voiceBits); 
extraBits = mod(-totalBits, bitspersymbol);
voiceBits = [voiceBits; zeros(extraBits,1)];



chunkedVoiceBits = reshape(voiceBits, bitspersymbol, [])'; %organizes  bits into a matrix (n bits per column)
voiceSymbols = BinarytoDecimal(chunkedVoiceBits); % Convert to symbols 
modVoice = exp(1j * 2*pi * (voiceSymbols) / M); % PSK Modulation

% Demodulation 
rxVoiceSymbols = zeros(length(modVoice), 1);
refPhases = exp(1j * 2 * pi * (0:M-1)/M);

for i = 1:length(modVoice)
    [~, val] = min(abs(modVoice(i) - refPhases));
    rxVoiceSymbols(i) = val - 1;
end
%----------------------------- 

%recovering voice message
rxVoiceBitsMat = rem(floor(rxVoiceSymbols(:) * ones(1, bitspersymbol) ./ 2.^(bitspersymbol-1:-1:0)),2); % converts decimal to binary
rxVoiceBits = rxVoiceBitsMat';  
rxVoiceBits = rxVoiceBits(:);
rxVoiceBits = rxVoiceBits(1:totalBits);    % remove zero padding
rxVoiceBitsBinMat = reshape(rxVoiceBits,8,[])';
rxVoiceDecimal = sum(rxVoiceBitsBinMat .* 2.^(7:-1:0), 2); % binary to decimal
rxVoice = (rxVoiceDecimal / 127.5) - 1;

disp('Playing Demodulated Voice Message...');
sound(rxVoice, fs);
pause(length(rxVoice)/fs + 1);

% Plotting Voice message
figure;
subplot(2,1,1);
plot(voice); title('Original Voice Signal');
subplot(2,1,2);
plot(rxVoice); title('Recovered Voice after PSK Modulation/Demodulation');

% PSK Constellation Plot (Voice message)
figure;
numPlot = min(length(modVoice), 5000); % limit points for display speed
plot(real(modVoice(1:numPlot)), imag(modVoice(1:numPlot)), '*');
axis equal;
grid on;
xlim([-1.2 1.2]); ylim([-1.2 1.2]);
title('PSK Constellation Diagram (Voice Message)');
xlabel('In-Phase'); ylabel('Quadrature');
