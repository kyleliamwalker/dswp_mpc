function [FFT] = fftofTS4(data,fs)
%FFT of Time Series
%data is the time series you want to analyse
%fs is the samping frequency

% Returns a structure, FFT, with:

    %.f - frequency bins [Hz]
    %.A - amplitude spectrum [m]
    %.A_mean - mean amplitude spectrum [m]
    %.S - energy spectrum for each gauge[m^2s]
    %.S_mean - mean energy spectrum [m^2s]

% fftofTS - J.Steynor original....

% S. Draycott edit 04/11/2016
   %- outputting S(f)
   %- removing requirement to be of length 2.
   % Making into structure
   % adding mean spectra
   
% Sorting size so FFT done over correct dimension
[a,b] = size(data);
if b > a
    data = data';
end
% FFT
% Length=max(size(data));
Length = 2^(nextpow2(max(size(data))));
nfft=Length;
DFT=fft(data,nfft)/(Length);
A=2*abs(DFT(1:nfft/2+1,:));
A(1)=A(1)./2;

% other outputs
f = (fs)./2*linspace(0,1,nfft/2+1)';
df=f(2)-f(1);
S=A.^2/(2*df);

%structure
FFT.S=S; FFT.A=A; FFT.f=f; FFT.DFT=DFT;
%means
FFT.S_mean=mean(S,2); FFT.A_mean=mean(A,2);

FFT.E = angle(DFT(1:nfft/2+1,:));
% FFT.phi = asin(DFT(1:nfft/2+1,:));
    
end
