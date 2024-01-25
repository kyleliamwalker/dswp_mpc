
function [ dist ] = ff_disturbances( robot, waves, jac, t, count )

sim_waves = waves;

[ sim_waves ] = evaluate_particles( t(count), robot.ekf.state(1), ...
            robot.ekf.state(2), sim_waves );
% [ sim_waves ] = evaluate_2nd_order_particles( t(count), robot.ekf.state(1), ...
%             robot.ekf.state(2), sim_waves );

[ sim_waves ] = pitch_particles( robot, sim_waves, t(count), jac );

v_b = jac' * [ sim_waves.vx; sim_waves.vz; sim_waves.v_theta ] ; 
a_b = jac' * [ sim_waves.ax; sim_waves.az; sim_waves.a_theta ];

M_A = robot.M_A;
D_q = robot.D_v_quad;
D_l = robot.D_v_lin;

% v_rob = [ robot.ekf.state(3:4); robot.state(6)];
% should this be relative or just particle velocity...?
dist = ( M_A*a_b + ( D_l + D_q.*abs(-v_b)) * (-v_b) );
% dist = ( M_A*awgn(a_b, 10) + ( D_l + D_q.*abs(-awgn(v_b,10))) * (-awgn(v_b,10)) );


end