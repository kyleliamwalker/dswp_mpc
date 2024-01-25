function [ robot ] = state_estimator( robot, tau, v_b, a_b, dt, count )

% need to incorporate jacobian somewhere

% measurement uses actual ground truth state
x = robot.state(1);     z = robot.state(2);
[ robot ] = ekf_2d(robot, [ x z ]', a_b, v_b, tau );

u_est = robot.ekf.x_hat(3);
w_est = robot.ekf.x_hat(4);

M_RB = robot.M_RB(1:2,1:2);
M_A = robot.M_A(1:2,1:2);
D_l = robot.D_v_lin(1:2,1:2);
D_q = robot.D_v_quad(1:2,1:2);

eqMotion = @(t, nu_in) -(M_RB + M_A) \ ...
    (M_A*a_b + ( D_l + D_q .* abs(nu_in-v_b) ) * (nu_in-v_b) ...
    - tau ) ;

nudot_est = eqMotion( [ 0 dt ], [ u_est, w_est ]' );

%% Plotting

robot.ekf.plot.x_hat(count, 1:4) = robot.ekf.x_hat;
robot.ekf.plot.y(count,:) = robot.ekf.y;

robot.ekf.plot.x_hat(count, 5) = nudot_est(1);
robot.ekf.plot.x_hat(count, 6) = nudot_est(2);

end