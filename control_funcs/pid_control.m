

function [ tau, gains, robot ] = pid_control( robot, jac, p_error, i_error, d_error, dt )

% [ p_error, i_error, d_error ] = calc_errors( robot, robot.p_error );
Kp = robot.Kp;      Ki = robot.Ki;      Kd = robot.Kd;

rot_p_error = jac' * p_error;
rot_i_error = jac' * i_error;
rot_d_error = jac' * d_error;

gains = Kp .* rot_p_error + Ki .* rot_i_error + Kd .* rot_d_error;
gains( gains >=  1) =  1;         gains( gains <= -1) = -1;
robot.mu = gains;

Tmax = robot.Tmax;                  %Max Thrust
tau_des = (1-exp(-dt/robot.t_m)) *  robot.tau_max .* gains;% + g_nu;        % quick fix?
mu = ucalloc(robot.K_u, robot.T_config, robot.T_weight, tau_des);
tau = Tmax * robot.T_config * mu;

% b_x = [ cos(robot.fA) cos(robot.fA) -cos(robot.aA) -cos(robot.aA) ];
% u_x = [ gains(1) gains(1) -gains(1) -gains(1) ];
% b_z = [ cos(robot.vA) cos(robot.vA) cos(robot.vA) cos(robot.vA) ];
% u_z = [ gains(2)+gains(3) gains(2)+gains(3) gains(2)-gains(3) gains(2)-gains(3) ];
% b_q = [ robot.l_x robot.l_x robot.l_x robot.l_x ];
% u_q = [ gains(2)+gains(3) gains(2)+gains(3) gains(2)-gains(3) gains(2)-gains(3) ];
% 
% tau = Tmax * [ b_x * u_x';
%                b_z * u_z';
%                b_q * u_q'];

end

