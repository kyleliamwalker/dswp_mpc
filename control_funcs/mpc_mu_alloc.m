

function [ tau, gains, robot ] = mpc_mu_alloc( robot, gains, dt )

gains( gains >=  1) =  1;         gains( gains <= -1) = -1;
robot.mu = gains;

Tmax = robot.Tmax;                  %Max Thrust
tau_des = (1-exp(-dt/robot.t_m)) *  robot.tau_max .* gains;% + g_nu;        % quick fix?
mu = ucalloc(robot.K_u, robot.T_config, robot.T_weight, tau_des);
tau = Tmax * robot.T_config * mu;

end

