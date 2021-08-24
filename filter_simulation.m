t = -1000:0.001:1000;
dt = t(2)-t(1);
f = linspace(-(1/(2*dt)),(1/(2*dt)),length(t));



%% Butterworth filter




%% Chebyshev filter




%% Elliptic(Cauer) filter




%% Bessel filter

% 10-th order
poly = 654729075 + 654729075*(1i*f) + 310134825*(1i*f).^2 + ...
    91891800*(1i*f).^3 + 18918900*(1i*f).^4 + 2837835*(1i*f).^5 + ...
    315315*(1i*f).^6 + 25740*(1i*f).^7 + 1485*(1i*f).^8 + 55*(1i*f).^9 + ...
    (1i*f).^10;
transfer = 654729075./poly;

% plot(f,log(abs(transfer)));
% plot(f,(transfer));

irf = ifft(ifftshift((transfer)));
plot(t,abs(irf));

%% Gaussian filter




%% Optimum L(Legender) filter




%% Linkwitz-riley filter



