% Author: Dr. Kyle L. Walker
% Description: Performs Deterministic Sea Wave Prediction for a given set
% of wave measurements and prediction parameters (time and distance). See
% papers listed in README.md for more information.

function [ t_p, spectra ] = dswp( TMeasure, x, w1  )

g = 9.81; % gravitational constant
t = linspace(0,TMeasure,length(w1)); %Measurement time
tf1 = max(t);
Ts = mean(diff(t)); %Time Interval between Measurements
m = length(w1); %Number of measurements taken
Fs = 1/Ts;  %Frequency of measurements
Fn = Fs/2;  %Niquist Frequency

N = pow2(nextpow2(m));  %Increase number of measurements to the next 
                        %largest power of 2 to increase efficiency of code

y = fft(w1,N); %calculate fourier transform of height measurements
Fv = linspace(0, 1, fix(N/2)+1)*Fn; %Remove conjugate Signal
Iv = 1:length(Fv);

%% Output Spectra

A = 1/N * 2.*abs(y(Iv))';
E = -angle(y(Iv))';

% filter components
A_check = A > 0.05*max(A);
f_check = Fv' > 0.05;
f_check_2 = Fv' < 0.5;

spectra.f = Fv(A_check & f_check & f_check_2)';
spectra.T = 1./spectra.f;

% Toggle whether components are noisy or not depending on desired test case
% Here for example, SNR = 15 

% spectra.A = A(A_check & f_check & f_check_2);
spectra.A = abs(awgn(A(A_check & f_check & f_check_2), 15, 'measured'));
spectra.H = 2.*spectra.A;

% spectra.E = E(A_check & f_check & f_check_2);
spectra.E = abs(awgn(E(A_check & f_check & f_check_2), 15, 'measured'));

% Calculate celerity of extreme frequency components
c_min = g/(2*pi*max(spectra.f));      
c_max = g/(2*pi*min(spectra.f));   %calculate min and max speed of wave
ts2 = x/c_min;
tf2 = tf1 + x/c_max;    %calculate prediction time region
if ts2 > tf2
    disp('Prediction Site beyond prediction diagram bounds');
   ts2 = tf1;
end

dt = t(2)-t(1);
t_p = (ts2:dt:tf2);
% output other wave parameters required
spectra.w = 2*pi*spectra.f;
spectra.k = spectra.w.^2 / g;
spectra.L = (2*pi)./spectra.k;

end
