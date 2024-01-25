

function [ robot ] = calc_errors_ekf( robot, prev_error, prev_ekf_error, count )

des = robot.ref_state;
state = robot.state;
ekf_state = [ robot.ekf.x_hat(1:2); robot.state(3); ...
              robot.ekf.x_hat(3:4); robot.state(6); ...
              robot.ekf.plot.x_hat(count, 5:6)'; robot.state(9) ];

robot.p_error = des - state(1:length(des));
robot.i_error = [ 0; 0; 0 ];
robot.d_error = robot.p_error - prev_error;

robot.ekf.p_error = des - ekf_state(1:length(des));
robot.ekf.i_error = [ 0; 0; 0 ];
robot.ekf.d_error = robot.ekf.p_error - prev_ekf_error;

% add i_error 

end

