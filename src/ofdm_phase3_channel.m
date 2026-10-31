
%% Phase 3 - OFDM Multipath Channel Estimation (Corrected)
clear; clc; close all;
rng(42);

%% Parameters
Nfft = 64;
Ncp = 16;
numFrames = 1000;
SNRdB = 0:4:24;

numBits = Nfft * numFrames * 2;

%% QPSK Data
txBits = randi([0 1], numBits, 1);
pairs = reshape(txBits, 2, []);

symbols = ((1-2*pairs(1,:)) + ...
    1j*(1-2*pairs(2,:))) / sqrt(2);

dataGrid = reshape(symbols, Nfft, numFrames);

%% Pilot and data OFDM symbols
pilot = ones(Nfft, 1);
X = zeros(Nfft, 2*numFrames);
X(:,1:2:end) = repmat(pilot, 1, numFrames);
X(:,2:2:end) = dataGrid;

%% OFDM modulation
txTime = ifft(X, Nfft, 1);
txCP = [txTime(end-Ncp+1:end,:); txTime];

%% Multipath channel
h = [0.9; 0.4+0.3j; 0.2-0.2j; 0; 0.1];
h = h / norm(h);

Htrue = fft([h; zeros(Nfft-length(h),1)], Nfft);

% Apply the channel to each CP-prefixed OFDM symbol.
% The CP is long enough to protect against this channel's memory.
rxChannel = complex(zeros(size(txCP)));

for n = 1:size(txCP,2)
    temp = conv(txCP(:,n), h);
    rxChannel(:,n) = temp(1:size(txCP,1));
end

% Calculate received signal power for noise setting
signalPower = mean(abs(rxChannel(:)).^2);

BER_noEQ = zeros(size(SNRdB));
BER_estEQ = zeros(size(SNRdB));
BER_perfectEQ = zeros(size(SNRdB));

%% SNR sweep
for k = 1:length(SNRdB)

    noisePower = signalPower / 10^(SNRdB(k)/10);
    noise = sqrt(noisePower/2) .* ...
        (randn(size(rxChannel)) + 1j*randn(size(rxChannel)));

    rxSignal = rxChannel + noise;

    %% Remove CP and FFT
    rxNoCP = rxSignal(Ncp+1:end,:);
    Y = fft(rxNoCP, Nfft, 1);

    Ypilot = Y(:,1:2:end);
    Ydata = Y(:,2:2:end);

    %% LS channel estimation from pilots
    Hest = Ypilot ./ repmat(pilot,1,numFrames);

    fprintf('SNR=%2d dB | Max channel error = %.6e\n', ...
        SNRdB(k), max(abs(mean(Hest,2)-Htrue)));

    %% Equalization
    ZnoEQ = Ydata;

    regularizer = 1e-8;
    Zest = Ydata .* conj(Hest) ./ ...
        (abs(Hest).^2 + regularizer);

    Hperfect = repmat(Htrue,1,numFrames);
    Zperfect = Ydata .* conj(Hperfect) ./ ...
        (abs(Hperfect).^2 + regularizer);

    %% BER
    detected = {ZnoEQ, Zest, Zperfect};
    results = zeros(1,3);

    for m = 1:3
        z = detected{m}(:);
        rxBits = zeros(2*length(z),1);
        rxBits(1:2:end) = real(z) < 0;
        rxBits(2:2:end) = imag(z) < 0;

        results(m) = mean(txBits ~= rxBits);
    end

    BER_noEQ(k) = results(1);
    BER_estEQ(k) = results(2);
    BER_perfectEQ(k) = results(3);

    fprintf(['No EQ=%.6f | LS EQ=%.6f | ' ...
        'Perfect EQ=%.6f\n\n'], results);
end

%% BER plot
figure;
semilogy(SNRdB, max(BER_noEQ,1/numBits), 'r-o', ...
    'LineWidth',1.5);
hold on;
semilogy(SNRdB, max(BER_estEQ,1/numBits), 'b-s', ...
    'LineWidth',1.5);
semilogy(SNRdB, max(BER_perfectEQ,1/numBits), 'g-^', ...
    'LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('OFDM over Multipath Channel - Corrected Model');
legend('Without EQ','LS Estimated EQ', ...
    'Perfect Channel EQ','Location','southwest');

%% Channel response plot
figure;
plot(0:Nfft-1, abs(Htrue), 'k-', 'LineWidth',1.5);
hold on;
plot(0:Nfft-1, abs(mean(Hest,2)), 'b--', 'LineWidth',1.2);
grid on;
xlabel('Subcarrier Index');
ylabel('|H(k)|');
title('Actual vs Average LS Channel Estimate');
legend('Actual Channel','Average LS Estimate');
