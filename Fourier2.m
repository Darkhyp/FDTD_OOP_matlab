f=10; %frequency of sine wave
overSampRate=30; %oversampling rate
fs=overSampRate*f; %sampling frequency
phase = 1/3*pi; %desired phase shift in radians
nCyl = 5; %to generate five cycles of sine wave
 
t=0:1/fs:nCyl*1/f; %time base
 
x=sin(2*pi*f*t+phase); %replace with cos if a cosine wave is desired
figure
    plot(t,x);
    title(['Sine Wave f=', num2str(f), 'Hz']);
    xlabel('Time(s)');
    ylabel('Amplitude');

%%
NFFT=1024; %NFFT-point DFT	 	 
X=fft(x,NFFT); %compute DFT using FFT	 	 
nVals=0:NFFT-1; %DFT Sample points	 	 
figure
    plot(nVals,abs(X));	 	 
    title('Double Sided FFT - without FFTShift');	 	 
    xlabel('Sample points (N-point DFT)')	 	 
    ylabel('DFT Values');

%%
NFFT=1024; %NFFT-point DFT	 	 
X=fft(x,NFFT); %compute DFT using FFT	 	 
nVals=(0:NFFT-1)/NFFT; %Normalized DFT Sample points	 	 
figure
    plot(nVals,abs(X));	 	 
    title('Double Sided FFT - without FFTShift');	 	 
    xlabel('Normalized Frequency')	 	 
    ylabel('DFT Values');

%%
NFFT=1024; %NFFT-point DFT	 	 
X=fftshift(fft(x,NFFT)); %compute DFT using FFT	 	 
fVals=(-NFFT/2:NFFT/2-1)/NFFT; %DFT Sample points	 	 
figure
    plot(fVals,abs(X));	 	 
    title('Double Sided FFT - with FFTShift');	 	 
    xlabel('Normalized Frequency')	 	 
    ylabel('DFT Values');

%%
NFFT=1024;	 	 
X=fftshift(fft(x,NFFT));	 	 
fVals=fs*(-NFFT/2:NFFT/2-1)/NFFT;	 	 
figure
    plot(fVals,abs(X),'b');	 	 
    title('Double Sided FFT - with FFTShift');	 	 
    xlabel('Frequency (Hz)')
    ylabel('|DFT Values|');

%%
NFFT=1024;
L=length(x);	 	 
X=fftshift(fft(x,NFFT));	 	 
Px=X.*conj(X)/(NFFT*L); %Power of each freq components	 	 
fVals=fs*(-NFFT/2:NFFT/2-1)/NFFT;	 	 
figure
    plot(fVals,Px,'b');	 	 
    title('Power Spectral Density');	 	 
    xlabel('Frequency (Hz)')	 	 
    ylabel('Power');
