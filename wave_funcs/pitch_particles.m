

function [ waves ] = pitch_particles( robot, waves, t, jac )

% x = robot.state(1);     z = robot.state(2);     theta = robot.state(3);
x = robot.ekf.x_hat(1);     z = robot.ekf.x_hat(2);     theta = robot.state(3);

strips = 10;
dx = robot.length/strips;
xvec = dx/2:dx:robot.length-dx/2 ; % body longitudinal location vector; centers
Lx = xvec - robot.length/2;

v_rot = zeros(3,numel(Lx));
a_rot = zeros(3,numel(Lx));  

v_theta = 0;     a_theta = 0;

for j = 1:length(Lx)

    x_local = x + cos(theta)*Lx(j);
    z_local = z + sin(theta)*Lx(j);

    [ loc_particles ] = evaluate_particles( t, x_local, ...
        z_local, waves ); 
%     [ loc_particles ] = evaluate_2nd_order_particles( t, x_local, ...
%           z_local, waves );

    v_rot(:,j) = jac' * [ loc_particles.vx; ...
        loc_particles.vz; 0 ];
    a_rot(:,j) = jac' * [ loc_particles.ax; ...
        loc_particles.az; 0 ];

    v_theta = v_theta + v_rot(2,j)/Lx(j);
    a_theta = a_theta + a_rot(2,j)/Lx(j);
    
end

waves.v_theta = v_theta;
waves.a_theta = a_theta;