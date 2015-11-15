%Matlab file used for class demo conducted during Session 7
%to illustrate the effects of aliasing AND the relationship
%between the DTFT and CTFT when the DT signal is obtained
%by sampling an analog CT signal.
%A word utterance is effectively sampled at four different
%rates and played back to hear the effects of aliasing.
%The DTFT of the signal obtained at each sampling rate is
%plotted.
clf
set(0,'defaultaxesfontsize',20);
%these commands read in the speech file
data=getspeech('enf1s1t0.wav');
Fs=12500;
%plot data to cut off silence
plot(data)
end1=input('end1?')
end2=input('end2?')
datar=data(end1:end2);
plot(datar)
[d,dsize]=size(datar);
%compute 8192 point FFT to view spectrum of speech signal
deltaf=Fs/8192;
freq=-Fs/2:deltaf:Fs/2-deltaf;
plot(freq,abs(fftshift(fft(datar,8192))),'linewidth',3)
xlabel('Analog Frequency (Hz)')
ylabel('Spectral Magnitude')
domega=2*pi/8192;
omega=-pi:domega:pi-domega;
pause
input('playback speech at original sampling rate of 12.5 KHz');
soundsc(datar,Fs)
input('reduce sampling rate by a factor of 2 thru decimation');
input('and playback speech at sampling rate of 6.25 KHz');
datar2=datar(1:2:dsize);
soundsc(datar2,Fs/2)
input('reduce sampling rate by a factor of 3 thru decimation');
input('and playback speech at sampling rate of 4.167 KHz')
datar3=datar(1:3:dsize);
soundsc(datar3,Fs/3)
input('reduce sampling rate by a factor of 4 thru decimation');
input('and playback speech at sampling rate of 3.125 KHz');
datar4=datar(1:4:dsize);
soundsc(datar4,Fs/4)
%compute via the FFT the DTFT of each of the
%four sampled versions of the original speech
%Plot magnitude of each DTFT over -pi to pi
subplot(221)
dr1=abs(fftshift(fft(datar,8192)));
plot(omega,dr1,'Linewidth',3)
axis([-pi pi 0 max(dr1)])
title('Fs=12.5 KHz')
xlabel('omega (radians)')
subplot(222)
dr2=abs(fftshift(fft(datar2,8192)));
plot(omega,dr2,'Linewidth',3)
axis([-pi pi 0 max(dr2)])
title('Fs=6.25 KHz')
xlabel('omega (radians)')
subplot(223)
dr3=abs(fftshift(fft(datar3,8192)));
plot(omega,dr3,'Linewidth',3)
axis([-pi pi 0 max(dr3)])
title('Fs=4.167 KHz')
xlabel('omega (radians)')
subplot(224)
dr4=abs(fftshift(fft(datar4,8192)));
plot(omega,dr4,'linewidth',3)
axis([-pi pi 0 max(dr4)])
title('Fs=3.125 KHz')
xlabel('omega (radians)')

