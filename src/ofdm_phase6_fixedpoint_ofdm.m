
%% PHASE 6: FIXED-POINT FFT IN AN OFDM SYSTEM
% Compare floating-point OFDM with fixed-point receiver FFT.
% Channel: AWGN
% Modulation: QPSK
% FFT: 64 points, CP: 16 samples

clear; clc; close all;
rng(10);

%% OFDM PARAMETERS
Nfft = 64;
Ncp = 16;
numSymbols = 2000;
SNRdB = 0:4:24;
wordLengths = [8 10 12 16];

fprintf('PHASE 6: FIXED-POINT FFT IN OFDM\n');
fprintf('================================\n');
fprintf('FFT length: %d | CP length: %d\n', Nfft, Ncp);
fprintf('OFDM symbols per SNR: %d\n\n', numSymbols);

%% GENERATE QPSK DATA
% Two bits per subcarrier: 1 -> positive, 0 -> negative.
txBits = randi([0 1], 2*Nfft*numSymbols, 1);

bitPairs = reshape(txBits, 2, []);
txQPSK = ((1 - 2*bitPairs(1,:)) + ...
    1j*(1 - 2*bitPairs(2,:))) / sqrt(2);

txQPSK = reshape(txQPSK, Nfft, numSymbols);

%% OFDM TRANSMITTER
txTime = ifft(txQPSK, Nfft, 1);

% Add cyclic prefix
txWithCP = [txTime(end-Ncp+1:end,:); txTime];

% Serialize OFDM frames
txSignal = txWithCP(:);

%% RESULT STORAGE
berFloat = zeros(size(SNRdB));
berFixed = zeros(length(wordLengths), length(SNRdB));

%% SNR SWEEP
for s = 1:length(SNRdB)

    % AWGN based on measured transmitted signal power
    signalPower = mean(abs(txSignal).^2);
    noisePower = signalPower / (10^(SNRdB(s)/10));

    noise = sqrt(noisePower/2) .* ...
        (randn(size(txSignal)) + 1j*randn(size(txSignal)));

    rxSignal = txSignal + noise;

    % Reshape received signal into OFDM frames
    rxFrames = reshape(rxSignal, Nfft+Ncp, numSymbols);

    % Remove cyclic prefix
    rxNoCP = rxFrames(Ncp+1:end,:);

    %% FLOATING-POINT RECEIVER
    rxFloat = fft(rxNoCP, Nfft, 1);
    rxFloatBits = qpskDemod(rxFloat(:));

    berFloat(s) = mean(rxFloatBits ~= txBits);

    %% FIXED-POINT RECEIVER
    for w = 1:length(wordLengths)
        B = wordLengths(w);

        rxFixed = complex(zeros(Nfft, numSymbols));

        for k = 1:numSymbols
            % Fixed-point FFT returns an unnormalized FFT.
            rxFixed(:,k) = fixedFFT6(rxNoCP(:,k), B);
        end

        rxFixedBits = qpskDemod(rxFixed(:));

        berFixed(w,s) = mean(rxFixedBits ~= txBits);
    end

    fprintf('SNR = %2d dB | Float BER = %.6g', ...
        SNRdB(s), berFloat(s));

    for w = 1:length(wordLengths)
        fprintf(' | %2d-bit BER = %.6g', ...
            wordLengths(w), berFixed(w,s));
    end
    fprintf('\n');
end

%% PLOT BER
figure;
semilogy(SNRdB, max(berFloat,1/(length(txBits))), ...
    '-ko', 'LineWidth', 1.8, 'MarkerSize', 6);
hold on;

plotColors = lines(length(wordLengths));

for w = 1:length(wordLengths)
    semilogy(SNRdB, ...
        max(berFixed(w,:),1/(length(txBits))), ...
        '-o', 'Color', plotColors(w,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6);
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('OFDM BER: Floating-Point vs Fixed-Point FFT');
legendStrings = [{'Floating-point'} ...
    arrayfun(@(b) sprintf('%d-bit fixed-point', b), ...
    wordLengths, 'UniformOutput', false)];
legend(legendStrings, 'Location', 'southwest');

%% PLOT FIXED-POINT BER ONLY
figure;
for w = 1:length(wordLengths)
    semilogy(SNRdB, ...
        max(berFixed(w,:),1/(length(txBits))), ...
        '-o', 'Color', plotColors(w,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6);
    hold on;
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('Effect of FFT Word Length on OFDM BER');
legend(arrayfun(@(b) sprintf('%d-bit', b), ...
    wordLengths, 'UniformOutput', false), ...
    'Location', 'southwest');

fprintf('\nPhase 6 completed.\n');

%% LOCAL FUNCTION: QPSK DEMODULATOR
function bits = qpskDemod(symbols)
    bits = zeros(2*length(symbols), 1);
    bits(1:2:end) = real(symbols) < 0;
    bits(2:2:end) = imag(symbols) < 0;
end

%% LOCAL FUNCTION: FIXED-POINT RADIX-2 FFT
function X = fixedFFT6(x, B)

    N = length(x);
    stages = log2(N);

    if 2^stages ~= N
        error('FFT length must be a power of two.');
    end

    % Signed fixed-point format:
    % Range approximately [-4, 4), with B-3 fractional bits.
    F = B - 3;

    % Quantize input
    x = qcomplex6(x, B, F);

    % Bit-reversal permutation
    bits = stages;
    xr = complex(zeros(N,1));

    for n = 0:N-1
        reversed = 0;

        for b = 1:bits
            reversed = 2*reversed + bitget(n,b);
        end

        xr(reversed+1) = x(n+1);
    end

    X = xr;

    % Radix-2 decimation-in-time FFT
    for stage = 1:stages

        M = 2^stage;
        halfM = M/2;

        for k = 1:M:N
            for j = 0:halfM-1

                twiddle = exp(-1j*2*pi*j/M);
                twiddle = qcomplex6(twiddle, B, F);

                upper = X(k+j);
                lower = X(k+j+halfM);

                product = qcomplex6(twiddle*lower, B, F);

                X(k+j) = qcomplex6(upper + product, B, F);
                X(k+j+halfM) = ...
                    qcomplex6(upper - product, B, F);
            end
        end
    end
end

%% LOCAL FUNCTION: QUANTIZE COMPLEX VALUES
function y = qcomplex6(x, B, F)

    step = 2^(-F);
    minValue = -2^(B-1)*step;
    maxValue = (2^(B-1)-1)*step;

    xr = min(max(real(x), minValue), maxValue);
    xi = min(max(imag(x), minValue), maxValue);

    y = round(xr/step)*step + ...
        1j*round(xi/step)*step;
end
