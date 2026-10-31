
%% PHASE 5: Radix-2 Fixed-Point FFT
% OFDM DSP Project
% Custom FFT with stage-wise quantization and saturation
% Base MATLAB only

clear;
clc;
close all;
rng(42);

%% 1. Parameters
N = 64;
wordLengths = [8 10 12 16];

fprintf('PHASE 5: FIXED-POINT RADIX-2 FFT\n');
fprintf('================================\n');
fprintf('FFT length: %d\n', N);
fprintf('FFT stages: %d\n\n', log2(N));

%% 2. Generate Complex Input
% Keep input components below full scale to avoid
% unnecessary clipping at the beginning.

x = randn(N,1) + 1j*randn(N,1);
x = 0.4*x/max([abs(real(x)); abs(imag(x))]);

% Quantize the input separately for each word length.

numFormats = length(wordLengths);

SQNR = zeros(numFormats,1);
EVM = zeros(numFormats,1);
maxError = zeros(numFormats,1);
overflowCount = zeros(numFormats,1);

%% 3. Run FFT for Each Fixed-Point Format
for idx = 1:numFormats

    B = wordLengths(idx);
    F = B-1;

    % Input quantization
    xq = qcomplex(x,B,F);

    % Our custom FFT scales by 1/2 at each stage.
    % Therefore compare against FFT(xq)/N.
    yFixed = fixedFFT(xq,B,F);

    yReference = fft(xq)/N;

    % Numerical error
    err = yReference-yFixed;

    signalPower = mean(abs(yReference).^2);
    errorPower = mean(abs(err).^2);

    SQNR(idx) = 10*log10(signalPower/errorPower);

    EVM(idx) = 100*sqrt( ...
        sum(abs(err).^2)/sum(abs(yReference).^2));

    maxError(idx) = max(abs(err));

    % Re-run while counting saturation events
    [~,overflowCount(idx)] = fixedFFT(xq,B,F);

    fprintf(['Word length: %2d bits | SQNR: %7.2f dB | ' ...
        'EVM: %8.4f %% | Max error: %.6g | ' ...
        'Saturations: %d\n'], ...
        B,SQNR(idx),EVM(idx),maxError(idx), ...
        overflowCount(idx));
end

%% 4. Plot SQNR
figure;
plot(wordLengths,SQNR,'o-','LineWidth',1.5);
grid on;
xlabel('Word Length (bits)');
ylabel('SQNR (dB)');
title('Fixed-Point FFT: SQNR vs Word Length');

%% 5. Plot EVM
figure;
plot(wordLengths,EVM,'s-','LineWidth',1.5);
grid on;
xlabel('Word Length (bits)');
ylabel('EVM (%)');
title('Fixed-Point FFT: EVM vs Word Length');

%% 6. Compare FFT Outputs at 8 and 16 Bits
y8 = fixedFFT(qcomplex(x,8,7),8,7);
y16 = fixedFFT(qcomplex(x,16,15),16,15);

figure;
plot(real(yReference),'k-o','LineWidth',1);
hold on;
plot(real(y8),'r--x','LineWidth',1);
plot(real(y16),'b-.','LineWidth',1);
grid on;
xlabel('FFT Bin');
ylabel('Real Amplitude');
title('Fixed-Point FFT vs Floating-Point Reference');
legend('Reference','8-bit FFT','16-bit FFT');

%% 7. Local Functions

function y = qcomplex(x,B,F)
    % Signed fixed-point quantization:
    % B total bits, F fractional bits.
    % Range: [-1, 1-2^(-F)] for the formats used here.

    step = 2^(-F);
    lo = -1;
    hi = 1-step;

    yr = round(real(x)/step)*step;
    yi = round(imag(x)/step)*step;

    yr = min(max(yr,lo),hi);
    yi = min(max(yi,lo),hi);

    y = complex(yr,yi);
end

function [y,overflowCount] = fixedFFT(x,B,F)

    N = length(x);
    stages = log2(N);

    if stages ~= floor(stages)
        error('FFT length must be a power of two.');
    end

    overflowCount = 0;

    %% Bit-reversal permutation
    bits = stages;
    xr = zeros(N,1);

    for n = 0:N-1
        reversed = 0;

        for b = 1:bits
            reversed = 2*reversed + ...
                bitget(n,b);
        end

        xr(reversed+1) = x(n+1);
    end

    a = qcomplex(xr,B,F);

    %% Radix-2 decimation-in-time FFT
    for s = 1:stages

        M = 2^s;
        halfM = M/2;

        for base = 1:M:N

            for j = 0:halfM-1

                idx1 = base+j;
                idx2 = base+j+halfM;

                % Twiddle factor W_M^j
                W = exp(-1j*2*pi*j/M);
                Wq = qcomplex(W,B,F);

                % Complex multiplication, quantized
                vRaw = a(idx2)*Wq;
                v = qcomplex(vRaw,B,F);

                u = a(idx1);

                % Divide by two at each stage to
                % control growth and avoid overflow.
                topRaw = (u+v)/2;
                botRaw = (u-v)/2;

                % Count component saturation events
                overflowCount = overflowCount + ...
                    sum(real(topRaw)<-1 | ...
                        real(topRaw)>1-2^(-F));
                overflowCount = overflowCount + ...
                    sum(imag(topRaw)<-1 | ...
                        imag(topRaw)>1-2^(-F));
                overflowCount = overflowCount + ...
                    sum(real(botRaw)<-1 | ...
                        real(botRaw)>1-2^(-F));
                overflowCount = overflowCount + ...
                    sum(imag(botRaw)<-1 | ...
                        imag(botRaw)>1-2^(-F));

                a(idx1) = qcomplex(topRaw,B,F);
                a(idx2) = qcomplex(botRaw,B,F);
            end
        end
    end

    y = a;
end
