clear sound;
[y,Fs]= audioread('eric.wav');   %y is the sampled data output and fs is the sample rate for this data

Y=fft(y); %fourier transform of the time domain signal "y" is "Y" in frequency domain
Yshifted=fftshift(Y); %this lets the signal Y to be centered at 0

%Plot the zero-frequency-shifted Fourier Transform "Yshifted"
n=length(Yshifted);
fshift = linspace( -Fs/2, Fs/2,length(Yshifted) );
figure(1);
subplot(3,1,1);
plot(fshift,Yshifted);
xlabel('Frequency');
title('Spectrum of audio file');

%Apply Ideal Filter -4k 4k
Yshifted(1: find(fshift>-4000,1))=0;
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

Am=max(abs(yTD)); %amplitude of the filtered msg that is centered at 0
Ac=2*Am;
fc=100000;

Fs_new=5*fc;
[Numer,Denom]=rat(Fs_new/Fs);
yTD=resample(yTD,Numer,Denom); 

t= linspace( 0,  length(yTD)/Fs_new,  length(yTD) );
fshift_new = linspace(-Fs_new/2, Fs_new/2,length(yTD));
t=transpose(t); %this is important to match yShiftedTD (transpose)

%DSB-SC modulation
s1=yTD.*cos(2*pi*fc*t); 

%Sketch the modulated DSB-SC in frequency domain
s1FD=real(fftshift(fft(s1)));
figure(2);
subplot(2,1,1);
plot(fshift_new,s1FD);
ylim([-3000 3000]);
xlabel('Frequency');
title('DSB-SC Frequency domain');

%DSB-TC modulation
s2=(Ac + yTD).*cos(2*pi*fc*t); 
%s2=ammod(real(yTD),fc,5*fc);

%Sketch the modulated DSB-TC in frequency domain
s2FD=real(fftshift(fft(s2)));
subplot(2,1,2);
plot(fshift_new,s2FD);
ylim([-3000 3000]);
xlabel('Frequency');
title('DSB-TC Frequency domain');

%Envelope detector for both modulation types 
envelope1=abs(hilbert(s1)); %%positive part of msg only
envelope2=abs(hilbert(s2)); %positive part of msg only

 %DC bias
TCdemodulated = detrend(envelope2); 

%Decrease the sampling frequency again, (Return to original fs)
envelope1 =resample(envelope1,Fs,5*fc);
TCdemodulated=resample(TCdemodulated,Fs,5*fc);

%Sound of envelope detectors
%sound(SCdemodulated,Fs);
%sound(real(TCdemodulated),Fs);

figure(3);
subplot(2,1,1);
plot(envelope1);
xlabel('Time');
title('Envelope Detection of DSB-SC ');

subplot(2,1,2);
plot(TCdemodulated);
xlabel('Time');
title('Envelope Detection of DSB-TC ');

%Detection of Am+SNR
out0 = awgn(s2,0);
envelope_out0=abs(hilbert(out0)); %positive part of msg only
 %DC bias
TCdemodulated_0 = detrend(envelope_out0); 
%Decrease the sampling frequency again, (Return to original fs)
TCdemodulated_0 =resample(TCdemodulated_0,Fs,5*fc);
%sound(real(TCdemodulated_0),Fs);
figure(4);
subplot(3,1,1);
plot(TCdemodulated_0);
xlabel('Time');
title('Envelope Detection of DSB-TC with SNR=0 ');

out10 = awgn(s2,10);
envelope_out10=abs(hilbert(out10)); %positive part of msg only
 %DC bias
TCdemodulated_10 = detrend(envelope_out10); 
%Decrease the sampling frequency again, (Return to original fs)
TCdemodulated_10 =resample(TCdemodulated_10,Fs,5*fc);
%sound(real(TCdemodulated_10),Fs);
figure(4);
subplot(3,1,2);
plot(TCdemodulated_10);
xlabel('Time');
title('Envelope Detection of DSB-TC with SNR=10 ');

out30 = awgn(s2,30);
envelope_out30=abs(hilbert(out30)); %positive part of msg only
 %DC bias
TCdemodulated_30 = detrend(envelope_out30); 
%Decrease the sampling frequency again, (Return to original fs)
TCdemodulated_30 =resample(TCdemodulated_30,Fs,5*fc);
%sound(real(TCdemodulated_30),Fs);
figure(4);
subplot(3,1,3);
plot(TCdemodulated_30);
xlabel('Time');
title('Envelope Detection of DSB-TC with SNR=30 ');

%coherent detection of dsb-sc +SNR
   % SNR=0 %
outS0 = awgn(s1,0);
coherentS0=outS0.*cos(2*pi*fc*t);
%Remove all frequencies greater than 4000 Hz
[b a] = butter (5, 4000.*2./Fs_new);
coherentFilS0 = filtfilt (b, a, coherentS0).*2; %Zero-phase digital filtering

% Decrease the sampling frequency again, (Return to original fs)
[Num,Den]=rat(Fs/Fs_new);
coherentFilS0 =resample(coherentFilS0,Num,Den);

%sound(real(coherentFilS0),Fs);

figure(5);
subplot(2,1,1);
plot(coherentFilS0);    %plot in time domain
xlabel('Time');
title('Snchronous detector of DSB-SC with SNR=0');
fshift = [0 fshift];
subplot(2,1,2);

plot(fshift,real(fftshift(fft(coherentFilS0))));    %plot in frequency domain
xlabel('Frequency');
title('Snchronous detector of DSB-SC with SNR=0');

    % SNR=10 %
outS10 = awgn(s1,10);
coherentS10=outS10.*cos(2*pi*fc*t);
%Remove all frequencies greater than 4000 Hz
coherentFilS10 = filtfilt (b, a, coherentS10).*2; %Zero-phase digital filtering

%Decrease the sampling frequency again, (Return to original fs)
[Num,Den]=rat(Fs/Fs_new);
coherentFilS10 =resample(coherentFilS10,Num,Den);

%sound(real(coherentFilS10),Fs);

figure(6);
subplot(2,1,1);
plot(coherentFilS10);    %plot in time domain
xlabel('Time');
title('Snchronous detector of DSB-SC with SNR=10');

subplot(2,1,2);
plot(fshift,real(fftshift(fft(coherentFilS10))));    %plot in frequency domain
xlabel('Frequency');
title('Snchronous detector of DSB-SC with SNR=10')

    % SNR=30 %
outS30 = awgn(s1,30);
coherentS30=outS30.*cos(2*pi*fc*t);
%Remove all frequencies greater than 4000 Hz
coherentFilS30 = filtfilt (b, a, coherentS30).*2; %Zero-phase digital filtering

%Decrease the sampling frequency again, (Return to original fs)
coherentFilS30 =resample(coherentFilS30,Num,Den);

%sound(real(coherentFilS30),Fs);

figure(7);
subplot(2,1,1);
plot(coherentFilS30);    %plot in time domain
xlabel('Time');
title('Snchronous detector of DSB-SC with SNR=30');

subplot(2,1,2);
plot(fshift,real(fftshift(fft(coherentFilS30))));    %plot in frequency domain
xlabel('Frequency');
title('Snchronous detector of DSB-SC with SNR=30');

% coherent detection with frequency error F=100.1 KHz instead of 100 KHz

fcerr=100100;
coherentErr=s1.*cos(2*pi*fcerr*t);
%Remove all frequencies greater than 4000 Hz
coherentErrFil = filtfilt (b, a, coherentErr).*2; %Zero-phase digital filtering

%Decrease the sampling frequency again, (Return to original fs)
coherentErrFil =resample(coherentErrFil,Num,Den);

%sound(real(coherentErrFil),Fs);

figure(8);
subplot(2,1,1);
plot(coherentErrFil);    %plot in time domain
xlabel('Time');
title('coherent detection with frequency error');

subplot(2,1,2);
plot(fshift,real(fftshift(fft(coherentErrFil))));    %plot in frequency domain
xlabel('Frequency');
title('Spectrum of coherent detection with frequency error')

%Repeat the coherent detection with phase error = 20 degree

Perr=(20*2*pi)/180;
coherentPErr=s1.*cos(2*pi*fc*t+Perr);
%Remove all frequencies greater than 4000 Hz
coherentPErrFil = filtfilt (b, a, coherentPErr).*2; %Zero-phase digital filtering

%Decrease the sampling frequency again, (Return to original fs)
coherentPErrFil =resample(coherentPErrFil,Num,Den);

%sound(real(coherentPErrFil),Fs);

figure(9);
subplot(2,1,1);
plot(coherentPErrFil);    %plot in time domain
xlabel('Time');
title('coherent detection with phase error 20');

subplot(2,1,2);
plot(fshift,real(fftshift(fft(coherentPErrFil))));    %plot in frequency domain
xlabel('Frequency');
title('Spectrum of coherent detection with phase error 20')
