function [ u ] = cas_opt_control_lin( cas_func, robot, v_b, a_b )


% sim_state = [ robot.state(4:6); robot.state(1:3) ];
sim_state = [ robot.ekf.x_hat(3:4); robot.state(6); ...
              robot.ekf.x_hat(1:2); robot.state(3)];

% noise included in ekf 
% state_noise = 0 * [ randn(1); randn(1); 0.1*randn(1); randn(1); randn(1); 0.1*randn(1) ];

x_ref = [ 0; 0; 0; robot.ref_state ];

% provides the solver with the starting point
p_vec(1:6,1) = sim_state;

for j = 1:robot.N+1
    
    p_vec(7:12,j) = [ v_b(:,j); a_b(:,j) ];
%     p_vec(13:18,j) = x_ref;
    p_vec(17:22,j) = x_ref;
    p_vec(13:16,j) = [ sim_state(6); sim_state(1:3) ];
    
end

u = full(cas_func(p_vec));


end

