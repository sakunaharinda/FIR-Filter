%H.L.S.H Jayasundara
%160243D

close all;
clear all;
%% Specifications
A = 2;
B = 4;
C = 3;

TAp = 0.05+(0.01*A); %Maximum passband ripple
TAa = 40+B; %Minimum stopband Attenuation
wp1 = (C*100)+400; %Lower passband edge
wp2 = (C*100)+700; %Upper passband edge
wa1 = (C*100)+250; %Lower stopband edge
wa2 = (C*100)+800; %Upper stopband edge
ws = 2*((C*100)+ 1200); %Sampling frequency

deltaP = (10^(0.05*TAp)-1)/(10^(0.05*TAp)+1); 
deltaA = 10^(-0.05*TAa);
delta = min(deltaP,deltaA);
Bt = min((wp1-wa1),(wa2-wp2));
wc1 = wp1 - Bt/2;  %Lower cutoff
wc2 = wp2 + Bt/2;  %Upper cutoff

%% Actual Attenuation
Aa = -20*log10(delta);
Ap = 20*log10((delta+1)/(1-delta));

%% Calculating Alpha
if Aa <= 21
    Alpha = 0;
elseif (21 < Aa)&& (Aa<=50)
    Alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
    Alpha = 0.1102*(Aa-8.7);
end

%% Calculating D
if Aa <= 21
    D = 0.9222;
else
    D = (Aa - 7.95)/14.36 ;
end

%% Finding N
N = ceil(((ws*D)/Bt) + 1);
if rem(N,2)==0
    N = N+1;
end

%% Impulse Response
nlim = (N-1)/2;
T = (2*pi)/ws;
n = (1 : nlim);
hn = (sin(wc2*T*n)- sin(wc1*T*n))./(n*pi);
hn_0 = 2*(wc2-wc1)/ws;
hn = cat(2,fliplr(hn),hn_0,hn);

figure;
stem((-nlim:nlim), hn,'fill');
title('Impulse Response Before going through the window h(nT)');
xlabel('Samples(n)');
ylabel('Amplitude');


%% Window Generation
n = -(N-1)/2:1:(N-1)/2; % length of the filter
beta = Alpha*sqrt(1-(2*n/(N-1)).^2);
lim = 100;
Ialpha = 1;
for k=1:lim
    temp = ((1./factorial(k))*(Alpha/2).^k).^2;
    Ialpha = Ialpha+temp;
end
Ibeta = 1;
for k=1:lim
    temp = ((1./factorial(k))*(beta/2).^k).^2;
    Ibeta = Ibeta+temp;
end
wkn = Ibeta/Ialpha;
figure;
stem(n,wkn,'fill');
title('Window function wk(nT)');
xlabel('Samples(n)');
ylabel('Amplitude');%% Magnitude response of the filter

%% Filter

hFilter = hn.*wkn;
[h,f] = freqz(hFilter);
fvtool(hFilter);  

%% Magnitude response of the filter in the passband

h = 20*log10(abs(h));
w = f/T;
figure
s = round(length(w)/(ws/2)*wc1);
f = round((length(w)/(ws/2)*wc2));
wpass = w(s:f);  %frequencies in passband
hpass = h(s:f);  %magnitude response in passband
plot(wpass,hpass)
axis([-inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Passband - Frequency Domain');

%% Filtering Signals
w1 = wc1/2;
w2 = (wc1+wc2)/2;
w3 = (ws/2 + wc2)/2;
n1 = 1:(1025 - N);
w = [w1 w2 w3];
xin = 0;
for w0=1:length(w)
   xin = xin + sin(w(w0).*n1.*T); %creating input signal
end
f1 = (1:1024)/1024*ws-ws/2;

Xjw = fft(xin,1024); %Input Response
Hjw = fft(hFilter,1024); %Filter response
Yjw = Xjw.*Hjw; %Filtered response
yout = ifft(Yjw); %Filtered Output
yo = sin(w2.*T.*(1:1024)); %Expected Output
Yojw = fft(yo,1024);
%% Time domain Representation

figure;
subplot(3,1,1) %Input Signal
stem(n1,xin);
axis([1 300 -3 3]);
title('Input signal (x(nT))');
xlabel('Samples (n)');
ylabel('Amplitude');

subplot(3,1,2)  %Output Signal
stem(f1,yout);
axis([1 300 -2 2]);
title('Output signal (y(nT))');
xlabel('Samples (n)');
ylabel('Amplitude');

subplot(3,1,3)  %Expected output signal
stem(f1,yo);
axis([1 300 -2 2]);
title('Expected Output signal (yo(nT))');
xlabel('Samples (n)');
ylabel('Amplitude');

%% Frequency domain Representation
figure;
Xk = fftshift(Xjw);
subplot(3,1,1);  %input signal
plot(f1,abs(Xk));
axis([-ws/2 ws/2 0 500]);
title('Input signal(X[k])');
ylabel('Amplitude');
xlabel('Frequency (rad/s)');

Yk = fftshift(Yjw);
subplot(3,1,2);  %output signal
plot(f1,abs(Yk));
title('Output signal(Y[k])');
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
axis([-ws/2 ws/2 0 500]);

Yok = fftshift(Yojw);
subplot(3,1,3);  %expected output signal
plot(f1,abs(Yok));
title('Expected Output signal(Yo[k])');
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
axis([-ws/2 ws/2 0 500]);
