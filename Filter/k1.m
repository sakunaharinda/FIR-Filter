m1 = 2;
m2 = 4;
m3 = 3;
%% Specifications

TAp = 0.05+0.01*m1; % dB max passband ripple
TAa = 40+m2; %dB min stopband attenuation
wp1 = m3*100+400; %rad/s lower passband edge
wp2 = m3*100+700; %rad/s upper passband edge
wa1 = m3*100+250; %rad/s lower stopband edge
wa2 = m3*100+800; %rad/s upper stopband edge
ws = 2*(m3*100+1200); %sampling freqency

bt1 = wp1-wa1; %lower transition width
bt2 = wa2-wp2; % upper transisiton width
bt = min(bt1,bt2); %critical transition width
wc1 = wp1-bt/2; % lower cutoff frequency
wc2 = wp2+bt/2; % upper cutoff frequency
T = 2*pi/ws; % sampling period

deltaP = (10^(0.05*TAp) - 1)/ (10^(0.05*TAp) + 1); % calculating delta
deltaA = 10^(-0.05*TAa);
delta = min(deltaP,deltaA);
Aa = -20*log10(delta); % Actual stopband attenuation

%% Calculating Alpha

if Aa<=21 % Calculating alpha
    alpha = 0;
elseif Aa>21 && Aa<= 50
    alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
    alpha = 0.1102*(Aa-8.7);
end

%% Calculating D

if Aa <= 21 % Calculating D
    D = 0.9222;
else
    D = (Aa-7.95)/14.36;
end

%% Choosing N

N = ceil(ws*D/bt +1); % order of the filter
if rem(N,2) == 0
    N = N+1;
end
n = -(N-1)/2:1:(N-1)/2; % length of the filter
beta = alpha*sqrt(1-(2*n/(N-1)).^2); %Calculating beta

%% Calculating I_alpha and I_beta

limit = 50;
Ialpha = 1;
for k = 1:limit
    temp = (1/factorial(k)*(alpha/2).^k).^2;
    Ialpha = Ialpha + temp;
end

Ibeta = 1;
k = 1;
for k = 1:limit
    temp = (1/factorial(k)*(beta/2).^k).^2;
    Ibeta = Ibeta + temp;
end

wknt = Ibeta/Ialpha;
figure
stem(n,wknt)
xlabel('n')
ylabel('Amplitude')
title('Kaiser Window - Time Domain');

%% Impulse Response Before window

nleft = -(N-1)/2:-1;
hntleft = 1./(nleft*pi).*(sin(wc2*nleft*T)-sin(wc1*nleft*T));
nright = 1:(N-1)/2;
hntright = 1./(nright*pi).*(sin(wc2*nright*T)-sin(wc1*nright*T));
hnt0 = 2/ws*(wc2-wc1);
hnt = [hntleft,hnt0,hntright];
figure
stem(n,hnt)
xlabel('n')
ylabel('Amplitude')
title(strcat(['Impulse Response (Before)']));

filter = hnt.*wknt;
fvtool(filter)
figure
[h,f] = freqz(filter);
w = f/T;
h = 20*log10(abs(h));
plot(w,h)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat(['Filter Response - Kaiser Window - Frequency Domain']));

figure
start = round(length(w)/(ws/2)*wc1);
finish = round((length(w)/(ws/2)*wc2));
wpass = w(start:finish);
hpass = (h(start:finish));
plot(wpass,hpass)
axis([-inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Passband - Frequency Domain');