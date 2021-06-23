clear all; clc; close all;

%% 1
fs = 10000;
N = 1024;
theta = 0.28808*pi;
n = 7;
[b1, a1] = butter(n, theta/pi);
[h1, w1] = freqz(b1, a1, N);

roots(b1); %zeros of lowpass
roots(a1); %polse of lowpass

%% 5
wp = 0.7*pi;
alfa = -(cos((theta + wp)/2)/cos((theta - wp)/2));

syms Z(z)
Z(z) = -(z.^-1 + alfa) / (1 + alfa*z.^-1);

%Nominator
pol = 0;
count = 0;
for i = b1
    pol = pol + i*Z(z).^count;
    count = count + 1;
end

%Denominator
den = 0;
count = 0;
for i = a1
    den = den + i*Z(z).^count;
    count = count + 1;
end

trfunc = pol/den;
[num, denum] = numden(trfunc);

bhp = sym2poly(num);
ahp = sym2poly(denum);
[hhp, whp] = freqz(bhp, ahp, N);

% Built-in 
[bbi, abi] = iirlp2hp(b1, a1, theta/pi, wp/pi);
[hbi, wbi] = freqz(bbi, abi, N);

figure
zplane(ahp,bhp);
title("calculated");

roots(bhp); %zeros of calculated transfer function
roots(ahp); %poles of calculated transfer function

figure
subplot(311); plot(w1/pi, 20*log10(abs(h1)));
title('LowPass From Qustion 1'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');
subplot(312); 
plot(wbi/pi, 20*log10(abs(hbi)));
title('Builtin HighPass'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');
subplot(313); 
plot(whp/pi, 20*log10(abs(hhp)));
title('Calculated HighPass'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');

%% 6
wp1 = 0.3*pi;
wp2 = 0.5*pi;

alfabp = cos((wp2 + wp1)/2) / cos((wp2 - wp1)/2);
kbp = cot((wp2 - wp1)/2) * tan(theta/2);


syms Z(z)
Zbp(z) = -(z.^-2 - (2*alfabp*kbp*z.^-1)/(kbp+1) + (kbp-1)/(kbp+1))/ (((kbp-1)/(kbp+1))*z.^-2 - (2*alfabp*kbp*z^-1)/(kbp+1) + 1);

%Denominator
polbp = 0;
count = 0;
for i = b1
    polbp = polbp + i*Zbp(z).^count;
    count = count + 1;
end

%Denominator
denbp = 0;
count = 0;
for i = a1
    denbp = denbp + i*Zbp(z).^count;
    count = count + 1;
end


trfuncbp = polbp/denbp;
[numbp, denumbp] = numden(trfuncbp);

bbp = sym2poly(numbp);
abp = sym2poly(denumbp);
[hbp, wbp] = freqz(bbp, abp, N);


% Built-in 
[abpbi, bbpbi] = iirlp2bp(b1, a1, theta/pi, [wp1/pi wp2/pi]);
[hbpbi, wbpbi] = freqz(abpbi, bbpbi, N);

figure
zplane(abp,bbp);
title("calculated");

roots(bbp); %zeros of calculated transfer function
roots(abp); %poles of calculated transfer function

figure
subplot(311); 
plot(w1/pi, 20*log10(abs(h1)));
title('LowPass From Qustion 1'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');
subplot(312); 
plot(wbpbi/pi, 20*log10(abs(hbpbi)));
title('Builtin BandPass'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');
subplot(313); 
plot(wbp/pi, 20*log10(abs(hbp)));
title('Calculated BandPass'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');

%% 7
A1 = 1;
A2 = 0.178; 
A3 = 10^-2; 

w1 = 0.4*pi;
w2 = 0.7*pi;

A = 40; 
beta = 3.395; 
Nmb = 2; 
M = 45; 
alfa = M/2;
G = [A1 A2];
w = [w1 w2];

h1 = zeros(M+1,1);
for i=1:M+1
    h1(i)=(G(1)-G(2))*(sin(w(1)*((i-1)-alfa))/(pi*((i-1)-alfa)))*(besselj(0,beta*(1-((i-1-alfa)/alfa)^2)^(1/2)))/besselj(0,beta); 
end

h2 = zeros(M,1);
for i=1:M+1
    h2(i)=(G(2)-A3)*(sin(w(2)*((i-1)-alfa))/(pi*((i-1)-alfa)))*(besselj(0,beta*(1-((i-1-alfa)/alfa)^2)^(1/2)))/besselj(0,beta);
end

hn = h1 + h2; 
[hfreq, wferq]=freqz(hn);

figure
stem(hn);
title('h[n]'); 
xlabel('n');

figure
plot(wferq/pi,20*log10(abs(hfreq)));
title('Magnitude of the frequency response of the designed DT-filter'); 
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (db)');

%% 8
n = 1:15;
n1 = 1:0.01:15;
x = cos(2*pi*0.3*n);
x1 = cos(2*pi*0.3*n1);
factor = 4;

pad = length(n)*(factor - 1);
half = ceil((length(x)+1)/2);
z = fft(x);

zeropading = zeros(1, pad);
zeropadded = [z(1:half) zeropading z(half+1:end)];

interpolated = real(ifft(zeropadded))*factor; 

figure
subplot(511);
plot(n1,x1);
axis([0 15 -inf inf]);
title("Original signal");
ylabel('Amplitude');
xlabel('Time');
subplot(512);
stem(n,x);
title("Sampled with n = 15");
ylabel('Amplitude');
xlabel('n');
subplot(513);
stem(abs(z));
title("Frequency Response");
ylabel('Amplitude');
xlabel('n');
subplot(514);
stem(abs(zeropadded));
title("Zero Padded Frequency Response with the factor 4");
ylabel('Amplitude');
xlabel('n');
subplot(515);
stem(interpolated);
title("Interpolated Signal Time Response");
ylabel('Amplitude');
xlabel('n');