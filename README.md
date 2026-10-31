# OFDM-Based Digital Wireless Communication System with Fixed-Point FFT Analysis

A MATLAB simulation of an OFDM wireless communication system that investigates BER performance under AWGN and multipath channel conditions, pilot-based channel estimation, one-tap equalization, and fixed-point FFT precision.

## Overview

Orthogonal Frequency Division Multiplexing (OFDM) is a multicarrier modulation technique used in modern wireless communication systems. This project implements a basic OFDM transceiver and progressively evaluates its performance under different channel conditions and numerical precision settings.

The project compares floating-point processing with a custom radix-2 fixed-point FFT software model to study how word length affects numerical accuracy and bit error rate (BER).

**Implementation:** MATLAB simulation  
**Modulation:** QPSK  
**FFT size:** 64 points  
**Cyclic prefix:** 16 samples  
**Channel models:** AWGN and a simulated multipath channel

## System Architecture

```text
Random Bit Generation
        |
        v
    QPSK Mapping
        |
        v
     OFDM IFFT
        |
        v
  Cyclic Prefix Addition
        |
        v
  AWGN / Multipath Channel
        |
        v
  Cyclic Prefix Removal
        |
        v
      FFT Receiver
        |
        v
 Pilot-Based Channel Estimation
        |
        v
   One-Tap Equalization
        |
        v
    QPSK Demodulation
        |
        v
      BER Analysis
```

The fixed-point FFT experiments model finite-precision arithmetic through quantization in a custom radix-2 FFT implementation.

## Project Modules

### 1. Basic OFDM Transceiver

- Generate random binary data and map it to QPSK symbols.
- Perform IFFT-based OFDM modulation.
- Add and remove the cyclic prefix.
- Recover the transmitted data using FFT-based demodulation.
- Evaluate bit recovery under an ideal channel.

### 2. BER Performance under AWGN

- Add complex additive white Gaussian noise at different SNR values.
- Measure BER across an SNR sweep.
- Plot BER to study the relationship between channel noise and receiver performance.

### 3. Multipath Channel and Channel Estimation

- Model a multipath channel using a finite impulse response.
- Insert known pilot OFDM symbols.
- Estimate the channel frequency response using the least-squares (LS) method.
- Apply one-tap frequency-domain equalization.
- Compare receivers without equalization, with LS channel estimation, and with perfect channel knowledge.

### 4. Fixed-Point Quantization

- Quantize OFDM time-domain samples at different word lengths.
- Calculate signal-to-quantization-noise ratio (SQNR) and error vector magnitude (EVM).
- Examine the effect of finite precision on the transmitted signal representation.

### 5. Custom Radix-2 Fixed-Point FFT

- Implement a 64-point radix-2 decimation-in-time FFT.
- Use bit-reversal input ordering and iterative butterfly stages.
- Quantize twiddle factors and intermediate results.
- Compare fixed-point FFT output with a floating-point MATLAB reference using SQNR, EVM, and maximum error.

### 6. Fixed-Point FFT in OFDM

- Integrate the custom FFT into the OFDM receiver.
- Compare BER for 8-, 10-, 12-, and 16-bit configurations.
- Evaluate the impact of numerical precision at different SNR values.

### 7. Final System Integration

- Combine OFDM modulation, cyclic prefix processing, multipath propagation, pilot-based LS channel estimation, equalization, and fixed-point FFT processing.
- Compare floating-point and fixed-point receiver BER.
- Visualize BER versus SNR and channel-estimation accuracy.

## Final Simulation Results

The following values are measured results from the final simulation. BER is the fraction of incorrectly detected bits among the simulated bits.

| SNR (dB) | Floating-point LS BER | Perfect-channel BER | 8-bit FFT BER | 10-bit FFT BER | 12-bit FFT BER | 16-bit FFT BER |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.289945 | 0.185289 | 0.293266 | 0.290352 | 0.289922 | 0.289930 |
| 4 | 0.169695 | 0.095398 | 0.175211 | 0.170195 | 0.169430 | 0.169617 |
| 8 | 0.080398 | 0.042789 | 0.086570 | 0.080523 | 0.080469 | 0.080438 |
| 12 | 0.037805 | 0.022117 | 0.043695 | 0.038313 | 0.037773 | 0.037844 |
| 16 | 0.020711 | 0.011914 | 0.026789 | 0.021039 | 0.020773 | 0.020734 |
| 20 | 0.010156 | 0.005180 | 0.015180 | 0.010953 | 0.010109 | 0.010148 |
| 24 | 0.004484 | 0.001789 | 0.008742 | 0.005039 | 0.004547 | 0.004484 |

**Interpretation:** BER decreases as SNR increases. The 8-bit configuration shows a larger BER penalty, while the 12- and 16-bit results are close to the floating-point LS baseline in this experiment. Perfect-channel equalization provides a reference using the known channel response rather than a pilot-derived estimate.

These values describe this particular simulation configuration; they are not universal performance guarantees.

## Tools and Concepts

- MATLAB and complex-valued signal processing
- OFDM modulation and demodulation
- QPSK modulation
- FFT and IFFT
- AWGN and multipath channel models
- Pilot-based least-squares channel estimation
- Frequency-domain equalization
- Fixed-point quantization
- Radix-2 FFT butterflies
- BER, SQNR, EVM and numerical error analysis

## How to Run

1. Clone or download this repository.
2. Open MATLAB or MATLAB Online.
3. Open the desired `.m` script.
4. Run the script and inspect the command-window results and generated figures.

The scripts use Base MATLAB functionality. Run the phases in order if you want to follow the development process.

## Limitations and Future Work

- The fixed-point FFT is a MATLAB software model, not a synthesized FPGA implementation.
- The current evaluation focuses on QPSK, a 64-point FFT and the channel models implemented in the scripts.
- The final BER results should be treated as simulation measurements, not as hardware benchmark results.
- Further work could include a bit-accurate RTL implementation, HDL simulation, FPGA resource and timing analysis, higher-order modulation, and more detailed channel models.

## Project Status

Completed at the MATLAB simulation level. FPGA RTL implementation and hardware validation have not been performed.

## Author

Manish Kumar K M Mani

GitHub: https://github.com/manishkumarkmmani-star  
LinkedIn: https://www.linkedin.com/in/manish-kumar-k-m-11897b379
