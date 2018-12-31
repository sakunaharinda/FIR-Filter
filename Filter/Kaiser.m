clear all;
close all;
%% Filter Specifications
A = 2;
B = 4;
C = 3;
TAp = 0.05+0.01*A; % Stopband attenuation.
TAa = 40+B; % Passband ripple.
wp1 = 400+C*100; % Lower passband frequency.
wp2 = 700+C*100; % Upper cutoff frequency.
wa1 = 250+C*100; % Lower stopband frequency.
wa2 = 800+C*100; % Upeer stopband frequency.
ws = 2*(1200+C*100); % Sampling frequency.
deltaP = (10^(0.05*TAp)-1)/(10^(0.05*TAp)+1);
deltaA = 10^(-0.05*TAa);
delta = min(deltaP,deltaA);
Bt = min((wp1-wa1),(wa2-wp2));
wc1 = wp1 - Bt/2;
wc2 = wp2 + Bt/2;
%% Actual Stopband Attenuation
Aa = -20*log10(delta); % Actual stopband ripple.
Ap = 20*log10((delta+1)/(1-delta)); % Actual Passband ripple.
%% Calculating Alpha
if Aa<=21
Alpha = 0;
elseif Aa <= 50
Alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
Alpha = 0.1102*(Aa-8.7);
end
%% Calculating D
if Aa <= 21
D = 0.9222;
else
D = (Aa - 7.95)/14.36 ;
end%% Calculating the Order(N)
N = ceil((ws*D/Bt)+1);
if rem(N,2)==0
N=N+1;
end
%% Calculating Impulse Response(h(n))
nMax = (N-1)/2;
T = (2*pi)/ws;
n = (1 : nMax);
hn = (sin(wc2*T*n)- sin(wc1*T*n))./(n*pi);
hn0 = 2*(wc2-wc1)/ws;
hn = cat(2,fliplr(hn),hn0,hn);
%% Calculating lo(alpha)
temp = 1;
loAlpha = 1;
k = 1;
while temp > 10^(-6)
temp = ((1/factorial(k))*(Alpha/2)^k)^2;
loAlpha = loAlpha+temp;
k = k+1;
end
%% Calculating lo(beeta)
loBeta = [];
for x = 1:nMax
Beta = Alpha*sqrt(1-(2*x/(N-1))^2);
temp = 1;
loTemp = 1;
k = 1;
while temp > 10^(-6)
temp = ((1/factorial(k))*(Beta/2)^k)^2;
loTemp = loTemp+temp;
k = k+1;
end
loBeta = cat(2,loBeta,loTemp);
end
loBeta = cat(2,fliplr(loBeta),loAlpha,loBeta);
%% Calculte Wkn and plot the function
wkn = loBeta/loAlpha;
figure;
stem((-nMax:1:nMax),wkn,'fill');
title('Window function wk(nT)');
xlabel('Samples(n)');
ylabel('Amplitude');%% Magnitude response of the filter
hFilter = hn.*wkn;
[h,f] = freqz(hFilter);
fvtool(hFilter);
figure;
plot(f/T, 20 * log10(abs(h)));
ylim([-100 4]);
xlim([0 1500]);
title('Magnitude response of the filter');
xlabel('Angular Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;
%% Magnitude response of the filter in the passband
figure;
plot(f/T, 20 * log10(abs(h)));
title('Magnitude response of the filter in the passband');
xlim([wc1 wc2]);
ylim([-TAp/2 TAp/2]);
xlabel('Angular Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;
%% Phase response of the filter
figure;
plot(f/T,unwrap(angle(h)));
title('Phase Response of the Filter');
xlabel('Frequency(rad/s)');
ylabel('Phase(rad)');
%% Impulse response before going through the window h(nT)
figure;
stem((-nMax:nMax), hn,'fill');
title('Impulse Response Before going through the window h(nT)');
xlabel('Samples(n)');
ylabel('Amplitude');
%% Impulse Response of the filter
figure;
stem((-nMax:nMax), hFilter,'fill');
title('Impulse Response of the filter hw(nT)');
xlabel('Samples(n)');
ylabel('Amplitude');
grid on;
%% Filter a Signal
w1 = wa1/2; % signal 1 (frequency 450 rad/s )
w2 = (wp1+wp2)/2; % signal 2 (frequency 1150 rad/s )
w3 = 1600; % signal 3 (frequency 1600 rad/s )
xn0 = @(x)(sin(w1*x*T)+sin(w2*T*x)+sin(w3*T*x));
n1 = 1:(1025-N); % we take 1025-N samples
xn = xn0(n1); % Input signal
Hk0 = fft(hFilter,1024); % Furier transform of the filter
Xk0 = fft(xn,1024); % Furier transform of the signal
Yk0 = Xk0.*Hk0; % Furier transform of the output
yFilter = ifft(Yk0); % Output signal of the filter
yOut = sin(w2*T*(1:1024)); % Expected output signal
%% Frequency domain representation of the input signal
figure;
f1 = (1:1024)/1024*ws-ws/2;
Xk = fftshift(Xk0);
subplot(3,1,1);
plot(f1,abs(Xk));
axis([-ws/2 ws/2 0 500]);
title('Input signal(X[k])');
ylabel('Amplitude');
xlabel('Frequency (rad/s)');
%% Frequency domain representation of the output signal
Yk = fftshift(Yk0);
subplot(3,1,2);
plot(f1,abs(Yk));
title('Output signal(Y[k])');
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
axis([-ws/2 ws/2 0 500]);
%% Expected output signal
Yk2 = fft(yOut,1024);
Yk2 = fftshift(Yk2);
subplot(3,1,3);
plot(f1,abs(Yk2));
title('Expected output signal(Y[k])');
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
axis([-ws/2 ws/2 0 500]);
%% Time domain representation of the input signal
figure;
subplot(3,1,1);
stem(n1,xn);
axis([1 300 -3 3]);
title('Input signal (x(nT))');
xlabel('Samples (n)');
ylabel('Amplitude');%% Time domain representation of the Output signal
subplot(3,1,2)
stem(f1,yFilter);
axis([1 300 -1.5 1.5]);
title('Output signal (y(nT))');
xlabel('Samples (n)');
ylabel('Amplitude');
%% Time domain representation of the expected output signal
subplot(3,1,3)
stem(f1,yOut);
axis([1 300 -1.5 1.5]);
title('Expected output signal (y(nT))');
xlabel('Samples (n)');
ylabel('Amplitude');