%Matlab file used for class demo conducted during Session 21
%illustrating efficacy of M-Channel PR Filter Bank
%based on a 2-channel QMF.  The 4 channel PR Bank is
%applied to a word utterance sampled at 12.5 KHz.
clf
clear all
set(0,'defaultaxesfontsize',20);
%Set up M-channel DFT filter bank:
M=4;
h0=[1 1]; h1=[1 -1];
h00=[1 0 1];   h10=[1 0 -1];
H(1,:)=conv(h0,h00);  G(1,:)=H(1,:);
H(2,:)=conv(h0,h10);  G(2,:)=-H(2,:);
H(3,:)=conv(h1,h00);  G(3,:)=-H(3,:);
H(4,:)=conv(h1,h10);  G(4,:)=H(4,:);
%
data=getspeech('0af1s1t0.wav');
Fs=12500;
clf
end1=500;
end2=12000;
%plot(data)
%end1=input('end1?');
%end2=input('end2?');
x=data(end1:end2);
%change sampling rate to 10 KHz
%x=resample(xr,4,5);
%Fs=(4/5)*Fs;
input('Utterance played back at 10 KHz sampling rate'); 
soundsc(x,Fs)
%create x_0 [n] through x_M-1 [n]
for m=1:M;
W(m,:)=conv(x,H(m,:));
X(m,:)=W(m,1:M:length(W(m,:)));
end
%create y_0[n] and y_M-1 [n]
for m=1:M;
Z(m,:)=zeros(1,M*length(X(m,:)));
Z(m,1:M:length(Z(m,:)))=X(m,:);
Y(m,:)=conv(Z(m,:),G(m,:));
end
y=zeros(1,length(Y(1,:)));
for m=1:M;
y=y+Y(m,:);
end
y=real(y);
input('QMF Output played back at Fs=10 KHz'); 
soundsc(y,Fs)
%plot and compare DTFT's of x[n] and y[n]
domega=2*pi/8192;
omega=-pi:domega:pi-domega;
yf1=abs(fftshift(fft(x,8192)));
yf2=abs(fftshift(fft(y,8192)));
subplot(211)
plot(omega,yf1,'Linewidth',3)
axis([-pi pi 0 max(yf1)])
title('DTFT of original utterance');
xlabel('omega (radians/s)');
subplot(212)
plot(omega,yf2,'Linewidth',3)
axis([-pi pi 0 max(yf2)])
title('DTFT of output of QMF');
xlabel('omega (radians/s)');
%plot filter responses h[n]
h0f=abs(fftshift(fft(H(1,:),512)));
h1f=abs(fftshift(fft(H(2,:),512)));
h2f=abs(fftshift(fft(H(3,:),512)));
h3f=abs(fftshift(fft(H(4,:),512)));
input('Plot DTFT of h0[n] thru h3[n]')
clf
domega=2*pi/512;
omega=-pi:domega:pi-domega;
plot(omega,h0f,'b','Linewidth',4)
axis([-pi pi 0 max(h0f)])
hold on
plot(omega,h1f,'r','Linewidth',4)
plot(omega,h2f,'g','Linewidth',4)
plot(omega,h3f,'c','Linewidth',4)
legend('H0(w)','H1(w)','H2(w)','H3(w)');
title('Frequency Response of h0[n] and h1[n]');
xlabel('omega (radians/s)');
hold off
pause
%plot and compare DTFT's of x0[n] and x1[n]
x0=X(1,:);  x1=X(2,:);
domega=2*pi/4096;
omega=-pi:domega:pi-domega;
xf1=abs(fftshift(fft(x0,4096)));
xf2=abs(fftshift(fft(x1,4096)));
subplot(211)
plot(omega,xf1,'Linewidth',3)
axis([-pi pi 0 max(xf1)])
title('DTFT of x0[n]');
xlabel('omega (radians/s)');
subplot(212)
plot(omega,xf2,'Linewidth',3)
axis([-pi pi 0 max(xf2)])
title('DTFT of x1[n]');
xlabel('omega (radians/s)');
