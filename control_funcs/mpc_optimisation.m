

function [ u_input ] = mpc_optimisation( t, dt, robot, waves, old_input, count )

max_iterations = 1000;
sim_robot = robot;
sim_waves = waves;
N = robot.N;

tol = 0.3;

sim_state_plot = zeros( N+1, 3 ); 
sim_control_plot = zeros( N, 3 );
cost_plot = zeros( N, 3 );
opt_cost_plot = zeros( N, 3 );

dJdu = zeros( N, 3 ) + tol/100;

input_weight = 0.01;

if all( isnan(old_input)) == 1
    
    sim_robot.state = robot.state;
    
    for n = 1:N
        
        k = count + n - 1;
        [ sim_waves ] = evaluate_particles( t(k), sim_robot.state(1), ...
        sim_robot.state(2), sim_waves );
        
        jac = [ cos(sim_robot.state(3))  sin(sim_robot.state(3))  0;
            -sin(sim_robot.state(3)) cos(sim_robot.state(3))  0;
            0           0       1];
        
        v_b = jac' * [ sim_waves.vx; sim_waves.vz; 0 ];
        a_b = jac' * [ sim_waves.ax; sim_waves.az; 0 ];
        
        [ tau, mu ] = pid_control( sim_robot, jac, sim_robot.p_error, sim_robot.i_error, sim_robot.d_error );
        sim_robot = dynamics( sim_robot, v_b, a_b, tau, jac, dt );

        [ sim_robot ] = calc_errors( sim_robot, sim_robot.p_error );
        
        sim_state_plot(n+1, :) = sim_robot.state(1:3);
        sim_control_plot(n,:) = mu;

        J = ( sim_robot.state(1:3) - sim_robot.ref_state ).^2 + ...
                    input_weight * mu.^2 ;
        cost_plot(n,:) = J;
        
    end

else
    
    sim_robot.state = robot.state;
    sim_control_plot(1:end-1, :) = old_input;

    for n = 1:N-1
        
        k = count + n - 1;
        [ sim_waves ] = evaluate_particles( t(k), sim_robot.state(1), ...
        sim_robot.state(2), sim_waves );
        
        jac = [ cos(sim_robot.state(3))  sin(sim_robot.state(3))  0;
            -sin(sim_robot.state(3)) cos(sim_robot.state(3))  0;
            0           0       1];
        
        v_b = jac' * [ sim_waves.vx; sim_waves.vz; 0 ];
        a_b = jac' * [ sim_waves.ax; sim_waves.az; 0 ];   
        
        [ tau, mu ] = mpc_mu_alloc( sim_robot, sim_control_plot(n,:)' );
        sim_robot = dynamics( sim_robot, v_b, a_b, tau, jac, dt );

        [ sim_robot ] = calc_errors( sim_robot, sim_robot.p_error );
        
        sim_state_plot(n+1, :) = sim_robot.state(1:3);

        J = ( sim_robot.state(1:3) - sim_robot.ref_state ).^2 + ...
                    input_weight * mu.^2 ;
        cost_plot(n,:) = J;

    end

    k = k + 1;
    
    [ sim_waves ] = evaluate_particles( t(k), sim_robot.state(1), ...
        sim_robot.state(2), sim_waves );
        
    jac = [ cos(sim_robot.state(3))  sin(sim_robot.state(3))  0;
        -sin(sim_robot.state(3)) cos(sim_robot.state(3))  0;
        0           0       1];

    v_b = jac' * [ sim_waves.vx; sim_waves.vz; 0 ];
    a_b = jac' * [ sim_waves.ax; sim_waves.az; 0 ];

    [ tau, mu ] = pid_control( sim_robot, jac, sim_robot.p_error, sim_robot.i_error, sim_robot.d_error );
    sim_robot = dynamics( sim_robot, v_b, a_b, tau, jac, dt );

    [ sim_robot ] = calc_errors( sim_robot, sim_robot.p_error );

    sim_state_plot(end, :) = sim_robot.state(1:3);
    sim_control_plot(end,:) = mu;

    J = ( sim_robot.state(1:3) - sim_robot.ref_state ).^2 + ...
                input_weight * mu.^2 ;
    cost_plot(end,:) = J;

end

%% OPTIMISATION STARTS HERE

sim_control_plot( sim_control_plot >  1) =  1;
sim_control_plot( sim_control_plot < -1) = -1;

for iter = 1:max_iterations

    sim_robot.state = robot.state;
    opt_state_plot = sim_state_plot;

    opt_control_plot = sim_control_plot - 0.01 * sign(dJdu);
    opt_control_plot(:,3) = sim_control_plot(:,3);

    opt_control_plot( opt_control_plot >  1) =  1;
    opt_control_plot( opt_control_plot < -1) = -1;

    for n = 1:N
        
        k = count + n - 1;
        [ sim_waves ] = evaluate_particles( t(k), sim_robot.state(1), ...
        sim_robot.state(2), sim_waves );
        
        jac = [ cos(sim_robot.state(3))  sin(sim_robot.state(3))  0;
            -sin(sim_robot.state(3)) cos(sim_robot.state(3))  0;
            0           0       1];
        
        v_b = jac' * [ sim_waves.vx; sim_waves.vz; 0 ];
        a_b = jac' * [ sim_waves.ax; sim_waves.az; 0 ];   
        
        [ tau, mu ] = mpc_mu_alloc( sim_robot, opt_control_plot(n,:)' );
        sim_robot = dynamics( sim_robot, v_b, a_b, tau, jac, dt );

        [ sim_robot ] = calc_errors( sim_robot, sim_robot.p_error );
        
        opt_state_plot(n+1, :) = sim_robot.state(1:3);

        J = ( sim_robot.state(1:3) - sim_robot.ref_state ).^2 + ...
                    input_weight * mu.^2 ;
        opt_cost_plot(n,:) = J;

    end

    dJdu = ( opt_cost_plot - cost_plot ) ./ ( opt_control_plot - sim_control_plot );
    
    tol_check_x = sum(abs(dJdu(:,1))) <= tol ;
    tol_check_z = sum(abs(dJdu(:,2))) <= tol ;
    tol_check_theta = sum(abs(dJdu(:,3))) <= tol ;

    if tol_check_x && tol_check_z
        fprintf('optimised \n');
        break
    else
        fprintf('looping \n');
    end
    
    sim_state_plot = opt_state_plot;
    sim_control_plot = opt_control_plot;
    cost_plot = opt_cost_plot;

end

u_input = opt_control_plot;
fprintf('timestep: %d \n', count);
fprintf('optimisation steps taken: %d \n', iter);

disp('---------------------------------------------------');


