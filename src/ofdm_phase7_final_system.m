
%% PHASE 7: FINAL OFDM SYSTEM INTEGRATION
% QPSK + OFDM + Multipath + AWGN + LS Channel Estimation
% + One-Tap Equalization + Fixed-Point Radix-2 FFT
%
% Base MATLAB only. No toolboxes required.

clear; clc; close all;
rng(42);

%% 1. PARAMETERS
Nfft = 64;
Ncp = 16;
numFrames = 1000;
SNRdB = 0:4:24;
wordLengths = [8 10 12 16];

fprintf('PHASE 7: FINAL OFDM SYSTEM\n');
fprintf('==========================\n');
fprintf('FFT size: %d | CP: %d | Frames: %d\n\n', ...
    Nfft, Ncp, numFrames);

%% 2. QPSK DATA GENERATION
numBits = 2 * Nfft * numFrames;

txBits = randi([0 1], numBits, 1);
pairs = reshape(txBits, 2, []);

txSymbols = ((1 - 2*pairs(1,:)) + ...
    1j*(1 - 2*pairs(2,:))) / sqrt(2);

dataGrid = reshape(txSymbols, Nfft, numFrames);

%% 3. PILOT AND DATA GRID
% Each frame contains one pilot followed by one data symbol.
pilot = ones(Nfft, 1);

X = complex(zeros(Nfft, 2*numFrames));
X(:,1:2:end) = repmat(pilot, 1, numFrames);
X(:,2:2:end) = dataGrid;

%% 4. OFDM MODULATION AND CYCLIC PREFIX
txTime = ifft(X, Nfft, 1);
txCP = [txTime(end-Ncp+1:end,:); txTime];

%% 5. MULTIPATH CHANNEL
h = [0.9; 0.4+0.3j; 0.2-0.2j; 0; 0.1];
h = h / norm(h);

Htrue = fft([h; zeros(Nfft-length(h),1)], Nfft);

% CP-aware channel processing:
% Extract the useful received samples after the CP interval.
rxChannel = complex(zeros(Nfft, 2*numFrames));

for n = 1:2*numFrames
    channelOutput = conv(txCP(:,n), h);
    rxChannel(:,n) = ...
        channelOutput(Ncp+1:Ncp+Nfft);
end

signalPower = mean(abs(rxChannel(:)).^2);

%% 6. RESULT ARRAYS
numSNR = length(SNRdB);
numWords = length(wordLengths);

berFloatLS = zeros(1, numSNR);
berFloatPerfect = zeros(1, numSNR);
berFixedLS = zeros(numWords, numSNR);

channelError = zeros(1, numSNR);

%% 7. SNR SWEEP
for s = 1:numSNR

    noisePower = signalPower / 10^(SNRdB(s)/10);

    noise = sqrt(noisePower/2) .* ...
        (randn(size(rxChannel)) + ...
         1j*randn(size(rxChannel)));

    rxSignal = rxChannel + noise;

    %% Floating-point FFT
    Yfloat = fft(rxSignal, Nfft, 1);

    YpilotFloat = Yfloat(:,1:2:end);
    YdataFloat = Yfloat(:,2:2:end);

    %% LS CHANNEL ESTIMATION
    HestFloat = YpilotFloat ./ repmat(pilot,1,numFrames);

    % Average pilot estimates to assess channel-estimation accuracy.
    Haverage = mean(HestFloat, 2);
    channelError(s) = max(abs(Haverage - Htrue));

    %% FLOATING-POINT EQUALIZATION
    regularizer = 1e-8;

    ZfloatLS = YdataFloat .* conj(HestFloat) ./ ...
        (abs(HestFloat).^2 + regularizer);

    Hperfect = repmat(Htrue,1,numFrames);

    ZfloatPerfect = YdataFloat .* conj(Hperfect) ./ ...
        (abs(Hperfect).^2 + regularizer);

    bitsFloatLS = qpskDemod7(ZfloatLS);
    bitsFloatPerfect = qpskDemod7(ZfloatPerfect);

    berFloatLS(s) = mean(bitsFloatLS ~= txBits);
    berFloatPerfect(s) = mean(bitsFloatPerfect ~= txBits);

    %% FIXED-POINT RECEIVER FFT
    for w = 1:numWords

        B = wordLengths(w);

        Yfixed = complex(zeros(Nfft, 2*numFrames));

        for n = 1:2*numFrames
            Yfixed(:,n) = fixedFFT7(rxSignal(:,n), B);
        end

        YpilotFixed = Yfixed(:,1:2:end);
        YdataFixed = Yfixed(:,2:2:end);

        % Estimate the channel using the fixed-point pilot FFT.
        HestFixed = YpilotFixed ./ repmat(pilot,1,numFrames);

        % One-tap regularized equalizer.
        Zfixed = YdataFixed .* conj(HestFixed) ./ ...
            (abs(HestFixed).^2 + regularizer);

        bitsFixed = qpskDemod7(Zfixed);
        berFixedLS(w,s) = mean(bitsFixed ~= txBits);
    end

    %% PRINT RESULTS
    fprintf('SNR = %2d dB | Float LS = %.6g | Float perfect = %.6g', ...
        SNRdB(s), berFloatLS(s), berFloatPerfect(s));

    for w = 1:numWords
        fprintf(' | %d-bit = %.6g', ...
            wordLengths(w), berFixedLS(w,s));
    end

    fprintf('\n');
end

%% 8. PLOT BER COMPARISON
% The floor is only for plotting measured zero-error points.
berFloor = 1/numBits;

figure;
semilogy(SNRdB, max(berFloatLS,berFloor), ...
    'k-o','LineWidth',1.7);
hold on;

semilogy(SNRdB, max(berFloatPerfect,berFloor), ...
    'k--^','LineWidth',1.7);

plotColors = lines(numWords);

for w = 1:numWords
    semilogy(SNRdB, max(berFixedLS(w,:),berFloor), ...
        '-o','Color',plotColors(w,:), ...
        'LineWidth',1.5);
end

grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate');
title('Final OFDM System: Floating-Point vs Fixed-Point FFT');

legendStrings = [{'Floating-point LS EQ', ...
                  'Floating-point perfect EQ'}];

for w = 1:numWords
    legendStrings{end+1} = ...
        sprintf('%d-bit fixed-point LS EQ', wordLengths(w));
end

legend(legendStrings,'Location','southwest');

%% 9. PLOT CHANNEL ESTIMATION ERROR
figure;
semilogy(SNRdB, max(channelError,1e-12), ...
    'b-o','LineWidth',1.7);
grid on;
xlabel('SNR (dB)');
ylabel('Maximum Error in Average Channel Estimate');
title('Channel Estimation Accuracy vs SNR');

%% 10. PLOT CHANNEL FREQUENCY RESPONSE
figure;
plot(0:Nfft-1,abs(Htrue),'k-','LineWidth',1.7);
grid on;
xlabel('Subcarrier Index');
ylabel('|H(k)|');
title('Multipath Channel Frequency Response');

fprintf('\nPhase 7 integration completed.\n');
fprintf('Bits tested per SNR: %d\n', numBits);
fprintf('FFT word lengths tested: ');
fprintf('%d ',wordLengths);
fprintf('\n');

%% LOCAL FUNCTION 1: QPSK DEMODULATION
function bits = qpskDemod7(symbols)

    symbols = symbols(:);
    bits = zeros(2*length(symbols),1);

    bits(1:2:end) = real(symbols) < 0;
    bits(2:2:end) = imag(symbols) < 0;
end

%% LOCAL FUNCTION 2: FIXED-POINT RADIX-2 FFT
function X = fixedFFT7(x,B)

    N = length(x);
    stages = log2(N);

    if 2^stages ~= N
        error('FFT length must be a power of two.');
    end

    % Signed fixed-point format with approximately [-4,4) range.
    F = B-3;

    x = quantizeComplex7(x,B,F);

    %% Bit-reversal permutation
    xr = complex(zeros(N,1));

    for n = 0:N-1
        reversed = 0;

        for b = 1:stages
            reversed = 2*reversed + bitget(n,b);
        end

        xr(reversed+1) = x(n+1);
    end

    X = xr;

    %% Iterative radix-2 DIT FFT
    for stage = 1:stages

        M = 2^stage;
        halfM = M/2;

        for k = 1:M:N
            for j = 0:halfM-1

                W = exp(-1j*2*pi*j/M);
                W = quantizeComplex7(W,B,F);

                upper = X(k+j);
                lower = X(k+j+halfM);

                product = quantizeComplex7(W*lower,B,F);

                X(k+j) = quantizeComplex7( ...
                    upper+product,B,F);

                X(k+j+halfM) = quantizeComplex7( ...
                    upper-product,B,F);
            end
        end
    end
end

%% LOCAL FUNCTION 3: COMPLEX FIXED-POINT QUANTIZER
function y = quantizeComplex7(x,B,F)

    step = 2^(-F);

    minValue = -2^(B-1)*step;
    maxValue = (2^(B-1)-1)*step;

    xr = min(max(real(x),minValue),maxValue);
    xi = min(max(imag(x),minValue),maxValue);

    y = round(xr/step)*step + ...
        1j*round(xi/step)*step;
end
