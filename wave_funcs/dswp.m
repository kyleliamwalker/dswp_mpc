function [ t_p, spectra ] = dswp( TMeasure, x, w1  )

% copied from ryan 

g = 9.81;
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
% min_index = find(Fv > 0.05, 1);
% max_index = find(Fv < 0.5, 1, 'last');
% index = find(A > 0.05*max(A));

A_check = A > 0.05*max(A);
f_check = Fv' > 0.05;
f_check_2 = Fv' < 0.5;

% spectra.f = Fv(index)';
spectra.f = Fv(A_check & f_check & f_check_2)';
spectra.T = 1./spectra.f;
% spectra.A = A(index);
% spectra.A = A(A_check & f_check & f_check_2);
spectra.A = abs(awgn(A(A_check & f_check & f_check_2), 15, 'measured'));
spectra.H = 2.*spectra.A;
% spectra.E = E(index);
% spectra.E = E(A_check & f_check & f_check_2);
spectra.E = abs(awgn(E(A_check & f_check & f_check_2), 15, 'measured'));

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

% for i = 1:length(spectra.T)
%     [ spectra.L(i), spectra.k(i), spectra.w(i) ] = disper( d, spectra.T(i) );
% end
% 
% spectra.L = spectra.L';
% spectra.k = spectra.k';
% spectra.w = spectra.w';

spectra.w = 2*pi*spectra.f;
spectra.k = spectra.w.^2 / g;
spectra.L = (2*pi)./spectra.k;

end
