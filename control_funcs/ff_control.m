
function [ tau, gains, robot ] = ff_control( robot, jac, p_error, i_error, d_error, pred_dist, dt )

u_ff = pred_dist ./ robot.tau_max;
vel_state = [ robot.ekf.x_hat(3:4); robot.state(6) ];

% [ p_error, i_error, d_error ] = calc_errors( robot, robot.p_error );
Kp = robot.Kp;      Ki = robot.Ki;      Kd = robot.Kd;

rot_p_error = jac' * p_error;
rot_i_error = jac' * i_error;
rot_d_error = jac' * d_error;

% limit the position gain to max speed of the ROV - set to 1m/s here
% although max BROV2 is 1.5m/s just to be sure it's achievable.
pos_gains = Kp .* rot_p_error + Ki .* rot_i_error + Kd .* rot_d_error;
pos_gains( pos_gains >=  1) =  1;         pos_gains( pos_gains <= -1) = -1;

gains = 1.2 .* ( pos_gains - vel_state ) + u_ff;
gains( gains >=  1) =  1;         gains( gains <= -1) = -1;

% Single loop
% gains = Kp .* rot_p_error + Ki .* rot_i_error + Kd .* rot_d_error + u_ff;

gains( gains >=  1) =  1;         gains( gains <= -1) = -1;
robot.mu = gains;

tau_des = (1-exp(-dt/robot.t_m)) *  robot.tau_max .* gains;        % quick fix?
mu = ucalloc(robot.K_u, robot.T_config, robot.T_weight, tau_des);
tau = robot.Tmax * robot.T_config * mu;

end

