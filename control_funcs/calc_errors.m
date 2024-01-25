

function [ robot ] = calc_errors( robot, prev_error )

des = robot.ref_state;
state = robot.state;

robot.p_error = des - state(1:length(des));
robot.i_error = [ 0; 0; 0 ];
robot.d_error = robot.p_error - prev_error;

% add i_error 

end

