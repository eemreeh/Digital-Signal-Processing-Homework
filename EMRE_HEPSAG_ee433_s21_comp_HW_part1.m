clear all; clc; close all;
% DSP Homework 1
Fs = 10000;
N = 1024;

% Butterworth Lowpass Filter
[b1,a1] = butter(7,0.28808);
[h1 w1] = freqz(b1,a1,N);

% Kaiser Window Method 
beta = 0;
n = 12+2;
Wc = 0.325;
hh = fir1(n,Wc,'low',kaiser(n+1,beta),'noscale');
[h2 w2] = freqz(hh,1,N);

% FIR Filter
[n1,fo,ao,w] = firpmord([Fs/2*0.25 Fs/2*0.4],[1 0],[0.1 0.1],Fs);
n1 = n1 + 1;                 %to make it same as  specifications
b2 = firpm(n1,fo,ao,w);
[h3 w3] = freqz(b2,1,1024);

% Amplitude of Filters
figure
subplot(311);
plot(w1/pi,20*log10(abs(h1)));
title('Butterworth Lowpass Filter AMPLITUDE');
ylabel('Amplitude');
xlabel('Normalized Frequency (x pi rad/sample)');
subplot(312);
plot(w2/pi,20*log10(abs(h2)));
title("Kaiser Window Method AMPLITUDE");
ylabel('Amplitude');
xlabel('Normalized Frequency (x pi rad/sample)');
subplot(313);
plot(w3/pi,20*log10(abs(h3)));
title("FIR Filter AMPLITUDE");
ylabel('Amplitude');
xlabel('Normalized Frequency (x pi rad/sample)');

% Phase of Filters
figure
subplot(311);
plot(w1/pi, 360/(2*pi)*unwrap(angle(h1)));
title('Butterworth Lowpass Filter PHASE');
ylabel('Phase (degrees)');
xlabel('Normalized Frequency (x pi rad/sample)');
subplot(312);
plot(w2/pi, 360/(2*pi)*unwrap(angle(h2)));
title("Kaiser Window Method PHASE");
ylabel('Phase (degrees)');
xlabel('Normalized Frequency (x pi rad/sample)');
subplot(313);
plot(w3/pi, 360/(2*pi)*unwrap(angle(h3)));
title("FIR Filter PHASE");
ylabel('Phase (degrees)');
xlabel('Normalized Frequency (x pi rad/sample)');

% Impulse Response of Filters
figure
subplot(311);
impz(b1,a1);
title("Butterworth Lowpass Filter Impulse Response");
ylabel('Amplitude');
xlabel('Samples (n)');
subplot(312);
impz(hh,1);
title("Kaiser Window Method Impulse Response");
ylabel('Amplitude');
xlabel('Samples (n)');
subplot(313);
impz(b2,1);
title("FIR Filter Impulse Response");
ylabel('Amplitude');
xlabel('Samples (n)');