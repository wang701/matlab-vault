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
plot(w, X_w, w, abs(Xr(1,:)));
axis([0 2*pi 0 max(X_w(1,:))]);
figure(2);
plot(w, X_w, w, abs(Xr(2,:)));
axis([0 2*pi 0 max(X_w(1,:))]);
% d) 	Both cases when N = 51, 101 cause a certain degree of aliasing since the
% 	original signal has infinite length. When N = 101, the aliasing effect
%	is greatly reduced.
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
figure(3);
plot(n1, x1, n1, x_hat1);
axis([0 20]);
% Original x[n]
n2 = 0:(N(2,:)-1);
x2 = a.^(abs(n2-D(2,:)));
% Construct x_hat[n]
k2 = 0:N(2,:)-1;
wk2 = 2*pi*k2/N(2,:);
X_wk2 = (1-a^2)./(1-2*a*cos(wk2)+a^2).*exp((-j)*D(2,:)*wk2);
x_hat2 = ifft(X_wk2,101);
% Construct x_a[n]
figure(4);
plot(n2, x2);
axis([0 100]);

% Problem 7.30
% a)
f1 = 1/128; f2 = 5/128; fc = 50/128;
n = 0:255;
x_org = cos(2*pi*f1*n) + cos(2*pi*f2*n);
x_c = cos(2*pi*fc*n);
x_am = x_org.*x_c;
figure(5);
subplot(3,1,1);
plot(n, x_org);
subplot(3,1,2);
plot(n,x_c);
subplot(3,1,3);
plot(n,x_am);
% b) and c) and d)
X_am_128 = fft(x_am(1:128),128);
X_am_128_t = fft(x_am(1:100),128);
X_am_256 = fft(x_am(1:180),256);
figure(6);
subplot(3,1,1);
plot(0:127,X_am_128);
subplot(3,1,2);
plot(0:127,X_am_128_t);
subplot(3,1,3);
plot(0:255,X_am_256);
figure(7);
subplot(3,1,1);
plot(abs(X_am_128));
subplot(3,1,2);
plot(abs(X_am_128_t));
subplot(3,1,3);
plot(abs(X_am_256));
% e) See attached derivations
