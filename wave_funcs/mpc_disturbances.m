
function [ v_b_plot, a_b_plot ] = mpc_disturbances( robot, waves, jac, t, count )

sim_waves = waves;

for j = 1:robot.N+1
    
    k = count + j - 1;
%     [ sim_waves ] = evaluate_particles( t(k), robot.state(1), ...
%         robot.state(2), sim_waves );
    [ sim_waves ] = evaluate_particles( t(k), robot.ekf.x_hat(1), ...
        robot.ekf.x_hat(2), sim_waves );
%     [ sim_waves ] = evaluate_2nd_order_particles( t(k), robot.ekf.x_hat(1), ...
%           robot.ekf.x_hat(2), sim_waves );
    [ sim_waves ] = pitch_particles( robot, sim_waves, t(k), jac );

    v_b = jac' * [ sim_waves.vx; sim_waves.vz; sim_waves.v_theta ];
    a_b = jac' * [ sim_waves.ax; sim_waves.az; sim_waves.a_theta ];
    
    v_b_plot(:,j) = v_b;
    a_b_plot(:,j) = a_b; 
    
end


end