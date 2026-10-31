
%% Phase 4 - Fixed-Point Quantization Analysis
% OFDM DSP Project
% Requires Base MATLAB only

clear;
clc;
close all;
rng(42);

%% STEP 1: OFDM Parameters
Nfft = 64;
Ncp = 16;
numOFDMSymbols = 2000;
bitsPerSymbol = 2;

numBits = Nfft * numOFDMSymbols * bitsPerSymbol;

% Test different fixed-point word lengths
wordLengths = [4 6 8 10 12 16];

%% STEP 2: Generate Random Bits
txBits = randi([0 1], numBits, 1);

%% STEP 3: QPSK Modulation
pairs = reshape(txBits, 2, []);

I = 1 - 2*pairs(1,:);
Q = 1 - 2*pairs(2,:);

qpskSymbols = (I + 1j*Q)/sqrt(2);

X = reshape(qpskSymbols, Nfft, numOFDMSymbols);

%% STEP 4: OFDM Modulation
txTime = ifft(X, Nfft, 1);

txCP = [txTime(end-Ncp+1:end,:); txTime];

txSignal = txCP(:);

%% STEP 5: Fixed-Point Quantization
% Normalize both I and Q using one common scale.
% Signed fixed-point range is approximately [-1, 1).

scale = max([abs(real(txSignal)); abs(imag(txSignal))]);

normalizedSignal = txSignal / scale;

SQNR = zeros(size(wordLengths));
EVM = zeros(size(wordLengths));
BER = zeros(size(wordLengths));

for k = 1:length(wordLengths)

    B = wordLengths(k);

    % Signed B-bit quantizer with B-1 fractional bits
    fractionalBits = B - 1;
    step = 2^(-fractionalBits);

    % Quantize real and imaginary components
    qReal = round(real(normalizedSignal)/step)*step;
    qImag = round(imag(normalizedSignal)/step)*step;

    % Saturate to signed fixed-point limits
    minVal = -1;
    maxVal = 1-step;

    qReal = min(max(qReal,minVal),maxVal);
    qImag = min(max(qImag,minVal),maxVal);

    quantizedSignal = scale*(qReal + 1j*qImag);

    %% Quantization Error
    errorSignal = txSignal - quantizedSignal;

    signalPower = mean(abs(txSignal).^2);
    errorPower = mean(abs(errorSignal).^2);

    SQNR(k) = 10*log10(signalPower/errorPower);

    %% OFDM Receiver
    rxFrame = reshape(quantizedSignal, ...
        Nfft+Ncp, numOFDMSymbols);

    rxTime = rxFrame(Ncp+1:end,:);

    Y = fft(rxTime,Nfft,1);

    rxSymbols = Y(:);

    %% Error Vector Magnitude
    referenceSymbols = X(:);

    EVM(k) = 100*sqrt( ...
        sum(abs(rxSymbols-referenceSymbols).^2) / ...
        sum(abs(referenceSymbols).^2));

    %% QPSK Demodulation
    rxI = real(rxSymbols) < 0;
    rxQ = imag(rxSymbols) < 0;

    rxPairs = [rxI.'; rxQ.'];
    rxBits = rxPairs(:);

    BER(k) = mean(txBits ~= rxBits);

    fprintf(['Word Length = %2d bits | ' ...
        'SQNR = %7.2f dB | ' ...
        'EVM = %7.3f %% | BER = %.6f\n'], ...
        B,SQNR(k),EVM(k),BER(k));

end

%% STEP 6: Plot SQNR
figure;
plot(wordLengths,SQNR,'o-','LineWidth',1.5);
grid on;
xlabel('Fixed-Point Word Length (bits)');
ylabel('SQNR (dB)');
title('Quantization SQNR vs Word Length');

%% STEP 7: Plot EVM
figure;
plot(wordLengths,EVM,'s-','LineWidth',1.5);
grid on;
xlabel('Fixed-Point Word Length (bits)');
ylabel('EVM (%)');
title('OFDM EVM vs Word Length');

%% STEP 8: Plot BER
figure;
semilogy(wordLengths,max(BER,1/numBits), ...
    'd-','LineWidth',1.5);
grid on;
xlabel('Fixed-Point Word Length (bits)');
ylabel('Measured BER (floor at 1/Nbits)');
title('BER vs Fixed-Point Word Length');
