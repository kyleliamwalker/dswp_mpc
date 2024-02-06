% Author: Dr. Kyle L. Walker
% Description: Calculate particle velocities and accelerations at a given
% location and according to a given JONSWAP spectra, according to a 2nd
% order model. Does not account for shallow/deep categorisation of waves.

function [ spectra ] = evaluate_2nd_order_particles( t, x, z, spectra )

H = spectra.H;
w = spectra.w;
T = spectra.T;
k = spectra.k;
L = spectra.L;
E = spectra.E;
d = spectra.d;

g = 9.81; 
vx = zeros(1,numel(t)); vz = zeros(1,numel(t));
ax = zeros(1,numel(t)); az = zeros(1,numel(t));

for i = 1:numel(T)   
    if d / L(i) > 0.5
        %deep
        c = sqrt(g*L(i)/(2*pi));
    
    elseif d / L(i) < 0.05
        %shallow
        c = sqrt(g*d);
    else 
        %intermediate
        c = sqrt( (g/k(i)) * tanh(2*pi*d/L(i)) );
    end
    
    % double check against 1st order, something not right
    vx = vx - g*H(i)/(2*c) * cosh(k(i)*(z+d))/cosh(k(i)*d) * cos(k(i)*x - w(i)*t + E(i)) + ...
        3/16 * c * k(i)^2 * H(i)^2 * cosh(2*k(i)*(z+d))/sinh(k(i)*d)^4 * cos(2*(k(i)*x - w(i)*t + E(i)));
    vz = vz + g*H(i)/(2*c) * sinh(k(i)*(z+d))/cosh(k(i)*d) * sin(k(i)*x - w(i)*t + E(i)) + ...
        3/16 * c * k(i)^2 * H(i)^2 * sinh(2*k(i)*(z+d))/sinh(k(i)*d)^4 * sin(2*(k(i)*x - w(i)*t + E(i)));
    ax = ax - 3/8 * c * k(i)^2 * H(i)^2 * w(i) * cosh(2*k(i)*(z+d))/sinh(k(i)*d)^4 * sin(2*(k(i)*x - w(i)*t + E(i))) - ...
        g*H(i)*w(i)/(2*c) * cosh(k(i)*(z+d))/cosh(k(i)*d) * sin(k(i)*x - w(i)*t + E(i));
    az = az + 3/8 * c * k(i)^2 * H(i)^2 * w(i) * sinh(2*k(i)*(z+d))/sinh(k(i)*d)^4 * cos(2*(k(i)*x - w(i)*t + E(i))) - ...
        g*H(i)*w(i)/(2*c) * sinh(k(i)*(z+d))/cosh(k(i)*d) * cos(k(i)*x - w(i)*t + E(i));

end

%vx = vx -0.01;
spectra.vx = vx; spectra.vz = vz;
spectra.ax = ax; spectra.az = az;

end