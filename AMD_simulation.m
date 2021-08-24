clear
clc

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

%% make exponential

disp('make exponential...');

decay = heaviside(t).*exp(-t/ideal_tau);
decay = decay/max(decay); %% 아마 scale을 1로 normalization 하기위해 있는것 같은데 없어도 1로 표준화가 된다?

%% choose filter

disp('1. Butterworth filter');
disp('2. Type I Chebyshev filter');
disp('3. Type II Chebyshev filter');
disp('4. Elliptic(Cauer) filter');
disp('5. Bessel filter');
disp('6. Gaussian filter');
disp('7. Optimum L(Legender) filter');
disp('8. Linkwitz-riley filter');
choose = input('Choose filter number(1~8): ')

irf=choose_filter(choose);

%% convolution with filter

signal = conv(decay,irf,'same');
signal = signal/max(signal);

irf = scale*irf;
signal = scale*signal;

%% digitize

disp('digitize...');
dt_digit = 1/digit_sampling_frequency;
t_digit = (t(1):dt_digit:t(end))*ns;

irf_digit = round(interp1(t,irf,t_digit)); %% interp1 : 1-D data interpolation. 
signal_digit = round(interp1(t,signal,t_digit));

%% interpolation

disp('sp line interpolation...');
dt_interp = 1/interp_sampling_frequency;
t_interp = (t_digit(1):dt_interp:t_digit(end))*ns;

irf_interp = interp1(t_digit,irf_digit,t_interp,'spline');
signal_interp = interp1(t_digit,signal_digit,t_interp,'spline');

irf_interp = irf_interp/max(irf_interp);
signal_interp = signal_interp/max(signal_interp);

figure(2);
plot(t_interp,irf_interp,t_interp,signal_interp);

%% find lifetime. 
%% 1. mean 값을 어떻게 할지 정하고 
%% 2. 정한 mean 값을 계산, 두 값을 빼면 그게 life time.

disp('find lifetime by integral...');

%% 적분을 이용한 평균 계산. 끝부분 ringing 일어난 부분 짤라서 중간만 적분하게 수정해야 한다.

T_e = sum(t_interp.*signal_interp)/sum(signal_interp); % 한 -2.5에서 10까지
T_e0 = sum(t_interp.*irf_interp)/sum(irf_interp); % 한 -7.5에서 7.5까지
tau = T_e - T_e0;
disp(tau);

