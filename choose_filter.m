function irf = choose_filter(choose)

%% parameter & constant
ns = 1;
GHz = 1;
MHz = 1e-3;

t = (-100:0.001:100)*ns;

dt = t(2)-t(1);
len = length(t);
f = linspace(-(1/(2*dt)),(1/(2*dt)),len);

ideal_tau = 3*ns;

% 8bit resolution
bit_num = 8; % bit
scale = 2^bit_num; % define dynamic range & cutoff decibel

digit_sampling_frequency = 125*MHz;
interp_sampling_frequency = 100*GHz;


%% 1. Butterworth filter. 해결
if choose == 1
N=input('Enter the value of the N: ');
f_cutoff=digit_sampling_frequency./2;
butter = (1+(f./f_cutoff).^2*N).^(-0.5);
figure(1);
plot(f,butter);
axis([0 0.09 -0.1 1.5]);
irf_butter = abs(ifftshift(ifft(butter)));
irf = irf_butter/max(irf_butter);



%% 2. Type I Chebyshev filter. 해결?
elseif choose == 2
N=input('Enter the value of the N: ');
f=2*f/digit_sampling_frequency;
Cn=cosh(N*acosh(f));
chebyPoly=Cn;
chebyshev1=(1+(0.3*chebyPoly).^2).^(-1/2);
figure(1);
plot(f,chebyshev1);
axis([0 0.09 -0.1 1.5]);
irf_chebyshev1=abs(ifftshift(ifft(chebyshev1)));
irf=irf_chebyshev1./max(irf_chebyshev1);



%% 3. Type II Chebyshev filter ...주파수 영역은 뜬다
elseif choose == 3
N=input('Enter the value of the N: ');
f_cutoff=digit_sampling_frequency./2;
[z2,p2,k2] = cheby2(N,30,2*pi*f_cutoff,'s');
%order, passband ripple, stopband attenuation, edge frequncy 순서로.
[b2,a2] = zp2tf(z2,p2,k2);
[h2,w2] = freqs(b2,a2,length(f));
w2=w2./(2*pi);
figure(1);
plot(w2,mag2db(abs(h2)));
axis([0 0.1 -50 5]);
grid;
xlabel('Frequency (GHz)');
ylabel('Attenuation (dB)');
legend('chebyshev2');
el=abs(ifftshift(ifft(h2)));
irf=el./max(el);



%% 4. Elliptic(Cauer) filter ...주파수 영역은 뜬다
elseif choose==4
N=input('Enter the value of the N: ');
f_cutoff=digit_sampling_frequency./2;
[ze,pe,ke] = ellip(N,3,30,2*pi*f_cutoff,'s'); 
%order, passband ripple, stopband attenuation, edge frequncy 순서로.
[be,ae] = zp2tf(ze,pe,ke);
[he,we] = freqs(be,ae,length(f));
we=we./(2*pi);
figure(1);
plot(we,mag2db(abs(he)));
axis([0 0.1 -40 5]);
grid;
xlabel('Frequency (GHz)');
ylabel('Attenuation (dB)');
legend('elliptic');
el=abs(ifftshift(ifft(he)));
irf=el./max(el);





%% 5. Bessel filter 몇몇 숫자에서 NaN이 뜬다. 9이상의 홀수에서(7까지는 NaN 안뜸.)
elseif choose == 5;
% 10-th order
%poly = 654729075 + 654729075*(1i*f) + 310134825*(1i*f).^2 + ...
%    91891800*(1i*f).^3 + 18918900*(1i*f).^4 + 2837835*(1i*f).^5 + ...
%    315315*(1i*f).^6 + 25740*(1i*f).^7 + 1485*(1i*f).^8 + 55*(1i*f).^9 + ...
%    (1i*f).^10;
%bessel = 654729075./poly;
% plot(f,log(abs(bessel)));
% plot(f,(bessel));
N=input('Enter the value of the N: ');
f_cutoff=digit_sampling_frequency./2;
f=f./f_cutoff;
theta=((2/pi)^(0.5))*((1i*f).^(N+0.5)).*exp(1i*f).*besselk(N+0.5,1i*f);
theta0=factorial(2*N)./((2.^N).*factorial(N));
transfer=theta0./abs(theta);
f0=f.*f_cutoff;
figure(1);
plot(f0,transfer);
axis([0 0.09 -0.1 1.5]);
irf_bessel = abs(ifftshift(ifft(transfer)));
irf = irf_bessel/max(irf_bessel);



%% 6. Gaussian filter
elseif choose == 6;
cutoff_frequency = digit_sampling_frequency/2;
sigma_f = cutoff_frequency*cutoff_frequency/(log(scale)); %% 표준편차 시그마 인듯.
filter_f = exp(-f.*f/sigma_f);
irf_g = abs(ifftshift(ifft(filter_f))); %% ifft = inverse fast fourier transform
%% 왜인지 모르겠는데 매트랩으로 fft나 ifft하면 도메인이 거꾸로 나온다. 그래서 ifftshift로 제대로 교정해줌.
irf = irf_g/max(irf_g);



%% 7. Optimum L(Legender) filter





%% 8. Linkwitz-riley filter

end