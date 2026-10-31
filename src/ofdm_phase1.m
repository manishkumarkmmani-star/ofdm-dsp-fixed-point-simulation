
%% OFDM Transmitter and Receiver - Phase 1
% Project: OFDM-Based Wireless Transceiver with DSP Analysis
% Stage: Basic transmission through an ideal channel
% Requires: Base MATLAB only

clear;
clc;
close all;

%% STEP 1: System Parameters
Nfft = 64;          % Number of OFDM subcarriers
Ncp = 16;           % Cyclic prefix length
numOFDMSymbols = 100;

bitsPerSymbol = 2;  % QPSK carries 2 bits per symbol
numBits = Nfft * numOFDMSymbols * bitsPerSymbol;

fprintf('OFDM Phase 1 Simulation\n');
fprintf('-----------------------\n');

%% STEP 2: Generate Random Input Bits
txBits = randi([0 1], numBits, 1);

%% STEP 3: QPSK Modulation
% Map each pair of bits to a complex symbol.
bitPairs = reshape(txBits, 2, []);
I = 1 - 2*bitPairs(1,:);
Q = 1 - 2*bitPairs(2,:);

qpskSymbols = (I + 1j*Q) / sqrt(2);

% Arrange symbols across subcarriers.
X = reshape(qpskSymbols, Nfft, numOFDMSymbols);

%% STEP 4: IFFT - Frequency Domain to Time Domain
txTime = ifft(X, Nfft, 1);

%% STEP 5: Add Cyclic Prefix
% Copy the last Ncp samples to the beginning of each symbol.
cyclicPrefix = txTime(end-Ncp+1:end, :);
txFrame = [cyclicPrefix; txTime];

% Serialize the transmitted samples.
txSignal = txFrame(:);

%% STEP 6: Ideal Wireless Channel
% No noise, fading, or distortion in Phase 1.
rxSignal = txSignal;

%% STEP 7: Receiver - Reconstruct OFDM Frames
rxFrame = reshape(rxSignal, Nfft + Ncp, numOFDMSymbols);

% Remove the cyclic prefix.
rxTime = rxFrame(Ncp+1:end, :);

%% STEP 8: FFT - Time Domain to Frequency Domain
Y = fft(rxTime, Nfft, 1);

%% STEP 9: QPSK Demodulation
rxSymbols = Y(:);

rxI = real(rxSymbols) < 0;
rxQ = imag(rxSymbols) < 0;

rxBitPairs = [rxI.'; rxQ.'];
rxBits = rxBitPairs(:);

%% STEP 10: Evaluate Performance
numErrors = sum(txBits ~= rxBits);
BER = numErrors / numBits;

fprintf('FFT size              : %d\n', Nfft);
fprintf('Cyclic prefix length  : %d samples\n', Ncp);
fprintf('Number of OFDM symbols: %d\n', numOFDMSymbols);
fprintf('Total transmitted bits: %d\n', numBits);
fprintf('Bit errors            : %d\n', numErrors);
fprintf('Bit Error Rate (BER)  : %.6g\n', BER);

if numErrors == 0
    fprintf('PASS: All bits recovered correctly!\n');
else
    fprintf('CHECK: Errors detected. Inspect the receiver.\n');
end

%% STEP 11: Plot Results

% A. Transmitted time-domain waveform
figure;
plot(real(txSignal(1:250)));
grid on;
xlabel('Sample Index');
ylabel('Amplitude');
title('OFDM Transmitted Signal - Real Part');

% B. QPSK constellation at transmitter
figure;
plot(real(qpskSymbols(1:1000)), ...
     imag(qpskSymbols(1:1000)), '.');
grid on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Transmitted QPSK Constellation');

% C. Received constellation
figure;
plot(real(rxSymbols(1:1000)), ...
     imag(rxSymbols(1:1000)), '.');
grid on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);
xlabel('In-Phase');
ylabel('Quadrature');
title('Received Constellation - Ideal Channel');

% D. Compare the first 40 transmitted and received bits
figure;
stairs(0:39, txBits(1:40), 'LineWidth', 1.5);
hold on;
stairs(0:39, rxBits(1:40), '--', 'LineWidth', 1.5);
grid on;
ylim([-0.2 1.2]);
xlabel('Bit Index');
ylabel('Bit Value');
legend('Transmitted Bits', 'Received Bits');
title('Transmitted vs Received Bits');
