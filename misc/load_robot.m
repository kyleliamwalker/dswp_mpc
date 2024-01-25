

function [ robot ] = load_robot( init_pose, des_pose, t )

% pos, vel, acc
robot.state = [ init_pose; 0; 0; 0; 0; 0; 0 ];
robot.ref_state = des_pose;
robot.mu = [ 0; 0; 0 ];

% state errors initialisation
robot.p_error = [ 0; 0; 0 ];
robot.i_error = [ 0; 0; 0 ];
robot.d_error = [ 0; 0; 0 ];
robot.ekf.p_error = [ 0; 0; 0 ];
robot.ekf.i_error = [ 0; 0; 0 ];
robot.ekf.d_error = [ 0; 0; 0 ];

% dimensions

robot.length = 0.457; 
robot.width = 0.575; 
robot.height = 0.254;
robot.fA = deg2rad(45);        %Forward Thruster Angle
robot.aA = deg2rad(45);        %Aft Thruster Angle
robot.vA = deg2rad(0);        %vertical Thruster Angle
robot.Tmax = 35;              %max thrust in each DoF (for allocation)
robot.l_x = 0.12;
robot.l_y = 0.218;
robot.l_perp = 0.16;         % perp distance of propeller to Com (45deg)

%% hydrostatics

robot.mDry = 11.5;      m = robot.mDry;
robot.Ixx = 0.207;      Ixx = robot.Ixx;
robot.Iyy = 0.253;      Iyy = robot.Iyy;
robot.Izz = 0.377;      Izz = robot.Izz;

robot.M_RB = [ m*eye(2) [0; 0];
                [0,0]   Iyy ];
            
robot.W = m*9.81;
robot.B = 114.8;

robot.xG = 0;           xG = robot.xG;
robot.zG = 0;           zG = robot.zG;
robot.xB = 0;           xB = robot.xB;
robot.zB = -0.028;      zB = robot.zB;
robot.r_g = [ xG, zG ];
robot.r_b = [ xB, zB ];
        
%% hydrodynamics 

% TRO code values
robot.Xu_dot = 10.565;      Xu_dot = robot.Xu_dot;      
robot.Mq_dot = 0.65;        Mq_dot = robot.Mq_dot;   
robot.Xq_dot = -0.67;        Xq_dot = robot.Xq_dot;
robot.Zw_dot = 18.68;       Zw_dot = robot.Zw_dot;
robot.Xw_dot = 0;           Xw_dot = robot.Xw_dot;
robot.Zq_dot = 0;           Zq_dot = robot.Zq_dot; 
% 
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9655050
% robot.Xu_dot = 6.36;      Xu_dot = robot.Xu_dot;      
% robot.Mq_dot = 0.135;        Mq_dot = robot.Mq_dot;   
% robot.Xq_dot = -0.67;        Xq_dot = robot.Xq_dot;
% robot.Zw_dot = 18.68;       Zw_dot = robot.Zw_dot;
% robot.Xw_dot = 0;           Xw_dot = robot.Xw_dot;
% robot.Zq_dot = 0;           Zq_dot = robot.Zq_dot; 

robot.Zu_dot = robot.Xw_dot;    Zu_dot = robot.Zu_dot;
robot.Mu_dot = robot.Xq_dot;    Mu_dot = robot.Mu_dot;
robot.Mw_dot = robot.Zq_dot;    Mw_dot = robot.Mw_dot;

robot.M_A = [ Xu_dot, Xw_dot, Xq_dot;
              Zu_dot, Zw_dot, Zq_dot;
              Mu_dot, Mw_dot, Mq_dot];
          
% robot.Xu_lin = 4.03;        Xu_lin = robot.Xu_lin; 
% robot.Zw_lin = 5.18;        Zw_lin = robot.Zw_lin;
% % not large enough? causes oscillations when using 0.07 for some reason.
% % robot.Mq_lin = 0.07;        Mq_lin = robot.Mq_lin;
% robot.Mq_lin = 0.2;        Mq_lin = robot.Mq_lin;
% 
% robot.Xu_quad = 18.18;        Xu_quad = robot.Xu_quad; 
% robot.Zw_quad = 36.99;        Zw_quad = robot.Zw_quad;
% %robot.Mq_quad = 1.55;        Mq_quad = robot.Mq_quad;
% robot.Mq_quad = 3;        Mq_quad = robot.Mq_quad;

% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9655050

robot.Xu_lin = 13.7;        Xu_lin = robot.Xu_lin; 
robot.Zw_lin = 33;        Zw_lin = robot.Zw_lin;
robot.Mq_lin = 0.8;        Mq_lin = robot.Mq_lin;
robot.Xu_quad = 141;        Xu_quad = robot.Xu_quad; 
robot.Zw_quad = 190;        Zw_quad = robot.Zw_quad;
robot.Mq_quad = 0.47;        Mq_quad = robot.Mq_quad;

robot.D_v_lin= eye(3) .* [ Xu_lin, Zw_lin, Mq_lin, ];
robot.D_v_quad = eye(3) .* [ Xu_quad, Zw_quad, Mq_quad ];

%% control
robot.Kp = [ 1.2; 1.2; 0.5 ];
robot.Ki = [ 0; 0; 0 ];
robot.Kd = [ 3; 3; 1.0 ];

robot.t_m = 0.1;

fA = robot.fA; aA = robot.aA; vA = robot.vA;
l_x = robot.l_x; 

% these need to be the extended versions, feel like they'll be huge though?
robot.T_weight = eye(8);
robot.K_u = eye(8) .* robot.Tmax;
robot.tau_max = [ 2*robot.Tmax*(cos(fA) + cos(aA)); ...
                  4*robot.Tmax*cos(vA); ...
                  4*robot.Tmax*robot.l_x ];
% 1-4 horizontal (1,2 front, 3,4 back), 5-8 vertical (5,6 front, 7,8 back)
robot.T_config = [ cos(fA) cos(fA) cos(aA) cos(aA) 0 0 0 0; ...
                   0 0 0 0 cos(vA) cos(vA) cos(vA) cos(vA); ...
                   0 0 0 0 -l_x -l_x l_x l_x ];

%% EKF Initialisation

robot.ekf.noise = 0.1;                % measurement noise
robot.ekf.Qd = blkdiag(5,5);                   % process noise variance?
robot.ekf.Rd = blkdiag(5,5);                  % measurment noise variance?
robot.ekf.h  = 0.1; 	  % sampling time

% Initialise Kalman Filter
robot.ekf.x_prd = [ robot.state(1:2); robot.state(4:5) ];
robot.ekf.P_prd = diag([1 1 1 1]);
robot.ekf.x_hat = robot.ekf.x_prd;

% robot.ekf.x = [0; 0];
robot.ekf.plot.x_hat = [ zeros(1,numel(t)) + robot.state(1);
                         zeros(1,numel(t)) + robot.state(2);
                         zeros(1,numel(t)) + robot.state(3);
                         zeros(1,numel(t)) + robot.state(4);
                         zeros(1,numel(t)) + robot.state(5);
                         zeros(1,numel(t)) + robot.state(6) ]';
% robot.ekf.robotPlots.y = [ zeros(1,numel(t)) + robot.state(1);
%                            zeros(1,numel(t)) + robot.state(2);
%                            zeros(1,numel(t)) + robot.state(3)]';
             
robot.ekf.state = [ robot.state(1:2); robot.state(4:5) ];
               
end

