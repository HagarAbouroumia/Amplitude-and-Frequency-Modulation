clc;
clear sound;
[y,Fs]= audioread('eric.wav');   %y is the sampled data output and fs is the sample rate for this data

Y=fft(y); %fourier transform of the time domain signal "y" is "Y" in frequency domain
Yshifted=fftshift(Y); %this lets the signal Y to be centered at 0

%Plot the zero-frequency-shifted Fourier Transform "Yshifted"
n=length(Yshifted);
fshift = linspace( -Fs/2, Fs/2,length(Yshifted) );
figure(1);
subplot(3,1,1);
plot(fshift,real(Yshifted) );
xlabel('Frequency');
title('Spectrum of audio file');

%Apply Ideal Filter -4k 4k
Yshifted(1: find(fshift>-4000,1) )=0;
Yshifted(find(fshift>4000,1):end)=0;

subplot(3,1,2);
plot(fshift,real(Yshifted));
xlabel('Frequency');
title('Audio Spectrum after filter '); 

%Obtain the filtered signal in time domain
yShiftedTD=ifftshift(Yshifted); %rearranges the zero-frequency-shifted Fourier Transfrom "Yshifted" back to the original transform output (ifftshift undoes the result of fftshift)
yTD=ifft(yShiftedTD); 

time = linspace(0,length(yTD)/Fs, length(yTD));
subplot(3,1,3);
plot(time,yTD);
xlabel('Time');
title('Audio in Time Domain after filter (the original msg i will want to recieve) '); 

%Sound the filtered audio signal
%sound(real(yTD),Fs);
fc=100000;

Fs_new=5*fc;
[Numer,Denom]=rat(Fs_new/Fs);
yTD=resample(yTD,Numer,Denom); 

t= linspace(0,length(yTD)/Fs_new,length(yTD));
fshift = linspace(-Fs_new/2,Fs_new/2,length(yTD));
fshift=transpose(fshift);
t=transpose(t); 
%this is important to match yShiftedTD (transpose)
int_m=cumsum(yTD)/Fs_new;
Am=max(abs(yTD));
beta = 0.01;
kf=(beta*4000)/Am;
NBFM_Mod=cos(2*pi*fc*t + kf*2*pi*int_m); %%  Beta=kf *Am/B.W = 200 *0.2027 / 4000 = 0.01 (Beta is smaller than 1 so it will be NBFM)

%Sketch the modulated NBFM frequency domain
figure(2)
s2FD=fft(NBFM_Mod);
subplot(2,1,1);
plot(fshift,real(fftshift(s2FD)));
ylim([-3000 3000]);
xlabel('Frequency');
title('NBFM Frequency domain');
output=diff(NBFM_Mod)* 450 ;
envelope=abs(hilbert(output));

%DC bias
 demodulated=detrend(envelope);
 

%Decrease the sampling frequency again, (Return to original fs)
[Num,Den]=rat(Fs/Fs_new);
demodulated=resample(demodulated,Num,Den);

%demodulated=fmdemod(NBFM_Mod,fc,Fs_new,500);

subplot(2,1,2);
plot(demodulated);
ylim([-0.2 0.2]);
xlabel('time');
title('demodulator NBFM in time domain');