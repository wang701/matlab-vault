% ECE 538 MATLAB HOMEWORK 3
% Problem 7.29
% a)
a = 0.8;
w = linspace(0,2*pi-2*pi/1024,1024);
X_w = (1-a^2)./(1-2*a*cos(w)+a^2);
% b)
N1 = 101;
D1 = (N1-1)/2;
w1 = linspace(0,2*pi-2*pi/N1,N1);
Xr_1 = zeros(1,1024);
for k = 0:N1-1
	wk = 2*pi*k/N1;
	X_wk = (1-a^2)./(1-2*a*cos(wk)+a^2)*exp((-j)*D1*wk);
	Xr_1 = Xr_1+X_wk*((sin(N1*(w-wk)/2))./(N1*sin((w-wk)/2)).*exp(-j*(N1-1)*(w-wk)/2));
end	
plot(w, abs(Xr_1));
axis([0 2*pi 0 max(Xr_1)])
%N = [21;101];
%D = (N-1)/2;
%Xr = zeros(2,1024);
%for k = 0:N-1
	%wk(1,:) = 2*pi*k/N(1,:);
	%wk(2,:) = 2*pi*k/N(2,:);
	%X_wk(1,:) = (1-a^2)./(1-2*a*cos(wk(1,:))+a^2)*exp((-j)*D(1,:)*wk(1,:));
	%X_wk(2,:) = (1-a^2)./(1-2*a*cos(wk(2,:))+a^2)*exp((-j)*D(2,:)*wk(2,:));
	%Xr(1,:) = Xr(1,:)+X_wk(1,:)*((sin(N(1,:)*(w-wk(1,:))/2))./(N(1,:)*...
	%sin((w-wk(1,:))/2)).*exp(-j*(N(1,:)-1)*(w-wk(1,:))/2));
	%Xr(2,:) = Xr(2,:)+X_wk(2,:)*((sin(N(2,:)*(w-wk(2,:))/2))./(N(2,:)*...
	%sin((w-wk(1,:))/2)).*exp(-j*(N(2,:)-1)*(w-wk(2,:))/2));
%end	
%plot(w, abs(Xr(2,:)));
%axis([0 2*pi 0 max(Xr(2,:))]);
