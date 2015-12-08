% ECE 538 MATLAB HOMEWORK 3
% Problem 7.29
% a)
close all;
a = 0.8;
w = linspace(0,2*pi-2*pi/1024,1024);
X_w = (1-a^2)./(1-2*a*cos(w)+a^2);
% b) and c) See attached plots
N = [21;101];
D = (N-1)/2;
Xr = zeros(2,1024);
% Get reconstructed X for N = 21
for k = 0:N(1,:)-1
	wk(1,:) = 2*pi*k/N(1,:);
	X_wk(1,:) = (1-a^2)./(1-2*a*cos(wk(1,:))+a^2)*exp((-j)*D(1,:)*wk(1,:));
	Xr(1,:) = Xr(1,:)+X_wk(1,:)*((sin(N(1,:)*(w-wk(1,:))/2))./(N(1,:)*...
	sin((w-wk(1,:))/2)).*exp(-j*(N(1,:)-1)*(w-wk(1,:))/2));
end	
% Get reconstructed X for N = 101
for k = 0:N(2,:)-1
	wk(2,:) = 2*pi*k/N(2,:);
	X_wk(2,:) = (1-a^2)./(1-2*a*cos(wk(2,:))+a^2)*exp((-j)*D(2,:)*wk(2,:));
	Xr(2,:) = Xr(2,:)+X_wk(2,:)*((sin(N(2,:)*(w-wk(2,:))/2))./(N(2,:)*...
	sin((w-wk(2,:))/2)).*exp(-j*(N(2,:)-1)*(w-wk(2,:))/2));
end	
figure(1);
plot(w, abs(X_w), w, abs(Xr(1,:)));
axis([0 2*pi 0 max(X_w(1,:))]);
legend('|X(\omega)|', '|X_r|');
xlabel('w');
title(...
'Fig 7.29-1 Magnitude of X(\omega) vs. reconstructed X_r, N = 21');
figure(2);
plot(w, abs(X_w), w, abs(Xr(2,:)));
axis([0 2*pi 0 max(X_w(1,:))]);
legend('|X(\omega)|','|X_r|');
xlabel('w');
title(...
'Fig 7.29-2 Magnitude of X(\omega) vs. reconstructed X_r, N = 101');
% d) 	Both cases when N = 51, 101 cause a certain degree of aliasing since the
% 	original signal has infinite length. When N = 101, the aliasing effect
%	is greatly reduced to the point that one could not observe.
% e)
% Original x[n]
n1 = 0:(N(1,:)-1);
x1 = a.^(abs(n1-D(1,:)));
% Construct x_hat[n]
k1 = 0:N(1,:)-1;
wk1 = 2*pi*k1/N(1,:);
X_wk1 = (1-a^2)./(1-2*a*cos(wk1)+a^2).*exp((-j)*D(1,:)*wk1);
x_hat1 = ifft(X_wk1,21);
% Construct x_a[n]
x_a1 = a.^(abs(n1+N(1,:)-D(1,:))) + a.^(abs(n1-D(1,:))) + ...
	a.^(abs(n1-N(1,:)-D(1,:)));
figure(3);
plot(n1, x1, n1, x_hat1, n1, x_a1);
legend('x_t[n]','x_{hat}[n]','x_a[n]');
title('Fig 7.29-3 x[n] vs. IFFT of X(\omega) vs. aliased x_a[n] when N = 21');
xlabel('n');
% Original x[n]
n2 = 0:(N(2,:)-1);
x2 = a.^(abs(n2-D(2,:)));
% Construct x_hat[n]
k2 = 0:N(2,:)-1;
wk2 = 2*pi*k2/N(2,:);
X_wk2 = (1-a^2)./(1-2*a*cos(wk2)+a^2).*exp((-j)*D(2,:)*wk2);
x_hat2 = ifft(X_wk2,101);
% Construct x_a[n]
x_a2 = a.^(abs(n2+N(2,:)-D(2,:))) + a.^(abs(n2-D(2,:))) + ...
	a.^(abs(n2-N(2,:)-D(2,:)));
figure(4);
plot(n2, x2, n2, x_hat2, n2, x_a2);
legend('x_t[n]','x_{hat}[n]','x_a[n]');
title('Fig 7.29-4 x[n] vs. IFFT of X(\omega) vs. aliased x_a[n] when N = 101');
xlabel('n');

% Problem 7.30
% a)
f1 = 1/128; f2 = 5/128; fc = 50/128;
n = 0:255;
x_org = cos(2*pi*f1*n) + cos(2*pi*f2*n);
x_c = cos(2*pi*fc*n);
x_am = x_org.*cos(2*pi*fc*n);
figure(5);
subplot(3,1,1);
plot(n, x_org);
xlabel('n');
%title('Fig 7.30-1','x[n]=cos(2\pif1n)+cos(2\pif2n)');
subplot(3,1,2);
plot(n,x_c);
xlabel('n');
title('x_c[n]=cos(2\pif_cn)');
subplot(3,1,3);
plot(n,x_am);
xlabel('n');
title('x[n]=x[n]*x_c[n]');
% b) and c) and d)
X_am_128 = fft(x_am(1:128),128);
X_am_128_t = fft(x_am(1:100),128);
X_am_256 = fft(x_am(1:180),256);
figure(6);
subplot(3,1,1);
plot(2*pi*(0:127)/128,X_am_128);
%title('Fig 7.30-2','128-pt DFT of x_am, 0<=n<=127');
xlabel('w');
subplot(3,1,2);
plot(2*pi*(0:127)/128,X_am_128_t);
title('128-pt DFT of x_am, 0<=n<=99');
xlabel('w');
subplot(3,1,3);
plot(2*pi*(0:255)/256,X_am_256);
title('256-pt DFT of x_am, 0<=n<=179');
xlabel('w');
figure(7);
subplot(3,1,1);
plot(2*pi*(0:127)/128,abs(X_am_128));
%title('Fig 7.30-3','|X_{am}|, 0<=n<=127');
xlabel('w');
subplot(3,1,2);
plot(2*pi*(0:127)/128,abs(X_am_128_t));
title('|X_{am}|, 0<=n<=99');
xlabel('w');
subplot(3,1,3);
plot(2*pi*(0:255)/256,abs(X_am_256));
title('|X_{am}|, 0<=n<=179');
xlabel('w');
% e) See attached derivations
