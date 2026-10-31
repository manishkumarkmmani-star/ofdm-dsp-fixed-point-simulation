
%% OFDM Phase 2 - AWGN Channel and BER Analysis
clear;
clc;
close all;

rng(42);  % Reproducible simulation

%% System Parameters
Nfft = 64;
Ncp = 16;
numOFDMSymbols = 2000;
bitsPerSymbol = 2;

numBits = Nfft * numOFDMSymbols * bitsPerSymbol;
SNRdB = 0:2:20;

%% Generate Random Bits
txBits = randi([0 1], numBits, 1);

%% QPSK Modulation
bitPairs = reshape(txBits, 2, []);

I = 1 - 2*bitPairs(1,:);
Q = 1 - 2*bitPairs(2,:);

qpskSymbols = (I + 1j*Q)/sqrt(2);

%% Arrange Symbols into OFDM Frames
X = reshape(qpskSymbols, Nfft, numOFDMSymbols);

%% IFFT
txTime = ifft(X, Nfft, 1);

%% Add Cyclic Prefix
txFrame = [txTime(end-Ncp+1:end,:); txTime];
txSignal = txFrame(:);

%% AWGN Channel and Receiver
BER = zeros(size(SNRdB));

signalPower = mean(abs(txSignal).^2);

for k = 1:length(SNRdB)

    snrLinear = 10^(SNRdB(k)/10);
    noisePower = signalPower / snrLinear;

    % Complex Gaussian noise
    noise = sqrt(noisePower/2) * ...
        (randn(size(txSignal)) + 1j*randn(size(txSignal)));

    rxSignal = txSignal + noise;

    %% Receiver
    rxFrame = reshape(rxSignal, Nfft+Ncp, numOFDMSymbols);

    % Remove cyclic prefix
    rxTime = rxFrame(Ncp+1:end,:);

    % FFT
    Y = fft(rxTime, Nfft, 1);

    rxSymbols = Y(:);

    %% QPSK Demodulation
    rxI = real(rxSymbols) < 0;
    rxQ = imag(rxSymbols) < 0;

    rxBitPairs = [rxI.'; rxQ.'];
    rxBits = rxBitPairs(:);

    %% Calculate BER
    errors = sum(txBits ~= rxBits);
    BER(k) = errors / numBits;

    fprintf('SNR = %2d dB | Errors = %6d | BER = %.6f\n', ...
        SNRdB(k), errors, BER(k));

    % Save received symbols at 4 dB and 16 dB
    if SNRdB(k) == 4
        rxSymbolsLow = rxSymbols;
    elseif SNRdB(k) == 16
        rxSymbolsHigh = rxSymbols;
    end
end

%% Plot 1 - BER vs SNR
figure;
semilogy(SNRdB, BER, 'o-', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('OFDM QPSK Performance over AWGN Channel');

%% Plot 2 - Constellation at Low SNR
figure;
plot(real(rxSymbolsLow(1:1000)), ...
     imag(rxSymbolsLow(1:1000)), '.');
grid on;
axis equal;
xlim([-2 2]);
ylim([-2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Received QPSK Constellation at 4 dB SNR');

%% Plot 3 - Constellation at High SNR
figure;
plot(real(rxSymbolsHigh(1:1000)), ...
     imag(rxSymbolsHigh(1:1000)), '.');
grid on;
axis equal;
xlim([-2 2]);
ylim([-2 2]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Received QPSK Constellation at 16 dB SNR');
