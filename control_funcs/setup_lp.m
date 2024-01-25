function [ cas_func ] = setup_lp( robot, dt )

import casadi.*

%% initialise symbolically

pq0 = MX.sym('pq0'); vx0 = MX.sym('vx0'); vz0 = MX.sym('vz0'); vq0 = MX.sym('vq0');
px = MX.sym('px');  pz = MX.sym('pz');  pq = MX.sym('pq');
vx = MX.sym('vx');  vz = MX.sym('vz');  vq = MX.sym('vq');
x = [vx; vz; vq; px; pz; pq];

ux = MX.sym('ux');  uz = MX.sym('uz');  uq = MX.sym('uq');
u = [ux; uz; uq];

vp_x = MX.sym('vp_x');  vp_z = MX.sym('vp_z');  vp_q = MX.sym('vp_q');
ap_x = MX.sym('ap_x');  ap_z = MX.sym('ap_z');  ap_q = MX.sym('ap_q');
d = [vp_x; vp_z; vp_q; ap_x; ap_z; ap_q; pq0; vx0; vz0; vq0];

% http://ijariie.com/AdminUploadPdf/ROV_MODELS__HYDRODYNAMICS_AND_ASSOCIATED_SENSORS_ijariie9466.pdfn
% https://discovery.ucl.ac.uk/id/eprint/10054515/1/ROVcontrol%20oe%202018.pdf
M = robot.M_RB;
M_A = robot.M_A;
D_l = robot.D_v_lin;
D_q = robot.D_v_quad;
B = robot.tau_max;
 
% seems to work for the surge and heave but makes the pitch blow up at 2 or
% 3 points? for the rest of the simulation it works fine.
% lin_Dq = D_q*(d(8:10)-d(1:3)).*abs(d(8:10)-d(1:3)) + ...
%     2*D_q*(d(8:10)-d(1:3)) .* ( (x(1:3)-d(1:3)) - (d(8:10)-d(1:3)) );

% does x here need to be deltax to hold?
% lin_q = sin(d(7)) + cos(d(7))*(x(6)-d(7));  
% G = [ 0; 0; (robot.r_g(2)*robot.W-robot.r_b(2)*robot.B) * lin_q ];   
G = [ 0; 0; (robot.r_g(2)*robot.W-robot.r_b(2)*robot.B) * x(6) ]; 

% is this right?
ode = [ -(M+M_A)^-1 * ( M_A*d(4:6) + ... 
        (D_l)*(x(1:3)-d(1:3)) ...
        + G - (1-exp(-dt/robot.t_m)) * B.*u );
        x(1:3)];

f = Function('f', {x,u,d}, {ode});

%% integrator options

% T = 10;         % Time Horizon
% N = 100;        % Number of control intervals

intg_options = struct;
intg_options.tf = dt;      % integrate over
intg_options.simplify = true;
intg_options.number_of_finite_elements = 4;

dae = struct;
dae.x = x;
dae.p = [u;d];
dae.ode = f(x,u,d);

intg = integrator('intg', 'rk', dae, intg_options);

res = intg('x0', x, 'p', [u;d]);
x_next = res.xf;

F = Function('F',{x,u,d},{x_next},{'x','u','d'},{'x_next'});

% sim = F.mapaccum(N);
        
%% opti stack setup

opti = casadi.Opti();

Nc = robot.N;

x = opti.variable(6,Nc+1);
u = opti.variable(3,Nc);
p = opti.parameter(22,Nc+1);    % state (6), dist (6), ref_state (6)

Q = 250;
R = 1;
P = 1;

% add intermediate control penalty
% opti.minimize(Q*sumsqr(x(4:5,:)-p(16:17,:)) + 2*Q*sumsqr(x(6,:)-p(18,:)) ...
%     + R*sumsqr(u) + ...
%     P*sumsqr(diff(u(1,:)))+P*sumsqr(diff(u(2,:)))+P*sumsqr(diff(u(3,:))) );

opti.minimize(Q*sumsqr(x(1:6,:)-p(17:end,:)) + ...
    R*sumsqr(u) + ...
    P*sumsqr(diff(u(1,:)))+P*sumsqr(diff(u(2,:)))+P*sumsqr(diff(u(3,:))) );

% opti.minimize(Q*sumsqr(x(:,end-1)-p(13:end,:)) + ...
%     + R*sumsqr(u) + ...
%     P*sumsqr(diff(u(1,:)))+P*sumsqr(diff(u(2,:)))+P*sumsqr(diff(u(3,:))) ...
%     + Q*sumsqr(x(:,end-1)-p(13:end,:)) );


for k = 1:Nc
    opti.subject_to(x(:,k+1) == F(x(:,k),u(:,k),p(7:16,k)));
    
    % including intermediate cost
%     opti.subject_to(u(:,k+1) == (1-exp(-dt/robot.t_m)) * u(:,k));
    
end

opti.subject_to(-1<=u<=1);
opti.subject_to(x(:,1)==p(1:6,1));

opti.solver('sqpmethod', struct('qpsol','osqp'));

% opti.set_value(p,[0;0;1;1;3;3;]);
% sol = opti.solve();

opts = struct;
opts.qpsol = 'qrqp';
opts.print_header = false;
opts.print_iteration = false;
opts.print_time = false;
opts.qpsol_options.print_iter = false;
opts.qpsol_options.print_header = false;
opts.qpsol_options.print_info = false;
opti.solver('sqpmethod',opts);

%% generate function to map from parameters to optimal control

cas_func = opti.to_function('cas_func',{p},{u},{'p'},{'u_opt'});

end

