% Author: Dr. Kyle L. Walker
% Description: Loads waves based on a first-order linear wave theory model,
% with data being sourced from the online repository of Cefas. Accounts for
% the dispersion relation.

function [ spectra ] = dswp_load_waves( t, x, dx )

spectra.d = 54;         % seabed depth
spectra.theta = 0;      % as planar waves, directionality is zero
spectra.rho = 1025;     % sea water density
df = 0.005;             % freq. bin width
g = 9.81;               % gravitational constant
 
%%%%% JONSWAP Parameters: Tp = 11.1, Hs = 3.24
spectra.f = [ 0.03, 0.04, 0.05, 0.06, 0.069, 0.075, 0.08, 0.085, 0.09, ...
    0.095, 0.1, 0.105, 0.111, 0.12, 0.13, 0.14, 0.153, 0.17, 0.19, 0.21, ...
    0.23, 0.25, 0.27, 0.29, 0.325, 0.39, 0.47 ];
spectra.T = [ 33.33, 25, 20, 16.67, 14.55, 13.33, 12.5, 11.76, 11.11, ...
    10.53, 10, 9.52, 8.99, 8.33, 7.69, 7.14, 6.56, 5.88, 5.26, 4.76, ...
    4.35, 4, 3.7, 3.45, 3.08, 2.56, 2.13 ];
spectra.S = [ 0.0041, 0.0303, 0.0226, 0.0473, 0.222, 0.4404, 2.3637, ...
    5.2224, 19.6683, 14.3635, 7.6429, 4.0457, 3.7947, 5.9308, 7.3271, ...
    5.0053, 2.81, 1.5941, 1.8596, 0.9285, 0.8465, 0.5248, 0.4293, 0.3484, ...
    0.3302, 0.0599, 0.0443 ];

spectra.A = sqrt(2.*df.*spectra.S);
spectra.H = 2*spectra.A;
rng(3)
spectra.E = 2*pi*rand(1,length(spectra.f));

d = spectra.d;  H = spectra.H;  T = spectra.T;  E = spectra.E;
[ spectra.w, spectra.L, spectra.k ] = dispersion( d, T, g );
w = spectra.w; k = spectra.k;

% calculate temporal wave height over the simulation at both the
% measurement location and the prediction location for cross-checking.
spectra.eta = zeros(1, numel(t));
spectra.pred = zeros(1, numel(t));

for i = 1:numel(T)   
    spectra.eta = spectra.eta + H(i) / 2 * cos(k(i)*(x-dx) - w(i)*t + E(i));
    spectra.pred = spectra.pred + H(i) / 2 * cos(k(i)*x - w(i)*t + E(i)); 
end

% if missing parameters in the definition above throw an error
if numel(spectra.T) == numel(spectra.H) && numel(spectra.T) == numel(spectra.E)
    return
else
    disp('BAD WAVES!!');
end

end