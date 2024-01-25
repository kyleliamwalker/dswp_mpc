function [ robot ] = ekf_2d(robot, state, a_p, v_p, tau)

% extended kalman filter estimation for rov from simulated noisy
% measurements
% state is [ px, pz, vx, vz ], expansion for pitch tbc

mDry = robot.mDry;
mAdx = robot.Xu_dot;
mAdz = robot.Zw_dot;

Xu_l = robot.Xu_lin;      Zw_l = robot.Zw_lin;   
Xu_q = robot.Xu_quad;      Zw_q = robot.Zw_quad;  

M_RB = blkdiag(mDry, mDry);
M_A = [ mAdx, 0;
        0, mAdz ];
D_l = blkdiag(Xu_l, Zw_l );
D_q = blkdiag(Xu_q, Zw_q );

% state is passed in as column vector

e = robot.ekf.noise;                % measurement noise
h = robot.ekf.h;

x_prd = robot.ekf.x_prd;
P_prd = robot.ekf.P_prd;

% initialization of Kalman filter
Qd = robot.ekf.Qd;
Rd = robot.ekf.Rd;

w = 0.2*randn(1);                     % process noise
% E = [ 0 0 e e ]';

%% Filter

% disturbance 
d_f = M_A*a_p + (D_l + D_q.*abs(v_p))*v_p;

% C is measurement matrix?
Cd = [1 0 0 0 
      0 1 0 0 ]; 

% noise added to measurement
y = state + e * w;
%y = awgn(state, 60, 'measured');
%y = state;

% KF gain
K = P_prd * Cd' * inv( Cd * P_prd * Cd' + Rd );

% corrector
IKC = eye(4) - K*Cd;
P_hat = IKC * P_prd * IKC' + K * Rd * K';
x_hat = x_prd + K * (y - Cd * x_prd);

% discrete-time extended KF-model
f_hat = [ x_hat(3:end)
          -(M_RB+M_A)^-1 * ( M_A*a_p + (D_l + D_q.*abs(x_hat(3:end)-v_p)) * (x_hat(3:end)-v_p) ...
          - tau) ];
%f_d   = x_hat + h * ( f_hat + [vx 0]') ;
f_d   = x_hat + h * f_hat;

% Predictor (k+1)  
% Ad = I + h * A and Ed = h * E
% where A = df/dx is linearized about x = x_hat

% don't know if this linearization is right?
Ad   = [ eye(2)   eye(2)*h
         zeros(2)   diag(-(M_RB+M_A)^-1 * (1 + h*2*D_l*(x_hat(3:end)-v_p) + M_A*a_p)) ];
%Ed = h * E;   
% Ed = h * (E + [d_f' 0 0]');
Ed = h * [ diag(d_f)
           diag([e e])];

% calc new x_prd and P_prd
x_prd = f_d;
P_prd = Ad * P_hat * Ad' + Ed * Qd * Ed';

%% save everything

robot.ekf.x_prd = x_prd;
robot.ekf.P_prd = P_prd;
% robot.ekf.x = x;
robot.ekf.x_hat = x_hat;
robot.ekf.y = y;


end