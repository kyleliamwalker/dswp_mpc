
% Author: Dr. Kyle L. Walker
% Description: Launchpad for dynamically simulating a PID controlled ROV in
% 3DOF (surge, heave, pitch) under the influence of wave disturbances. #
% Desired set-point is altered by modifying ref_pose. An Extended Kalman
% Filter is employed for estimating the vehicle state. A GIF is produced at 
% the end of the simulation to visualise the vehicle motion.

%%
clear variables
%close all
clc

% set time parameters 
dt = 0.1;
t = 300:dt:500;

% set initial and desired rov pose
robot_pose = [ 50; -5; 0 ];
ref_pose = [ 50; -5; 0 ];

%% load wave parameters and vehicle specification
waves = load_waves( t, robot_pose(1), 1 );
% waves = load_2nd_order_waves( t, robot_pose(1) );

brov2 = load_robot( robot_pose, ref_pose, t );

%% run dynamic simulation
for k = 1:length(t)

    % calculate particle motions at each timestep
    [ waves ] = evaluate_particles( t(k), brov2.state(1), ...
        brov2.state(2), waves );
%     [ waves ] = evaluate_2nd_order_particles( t(k), brov2.state(1), ...
%             brov2.state(2), waves );
    waves.plot.vx(k) = waves.vx;    waves.plot.ax(k) = waves.ax;
    waves.plot.vz(k) = waves.vz;    waves.plot.az(k) = waves.az;

    % calculate jacobian
    jac = [ cos(brov2.state(3))  sin(brov2.state(3))  0;
            -sin(brov2.state(3)) cos(brov2.state(3))  0;
            0           0       1];
  
    % estimate pitch particles (alternative method)
    [ waves ] = pitch_particles( brov2, waves, t(k), jac );
    waves.plot.v_theta(k) = waves.v_theta;    waves.plot.a_theta(k) = waves.a_theta;
    
    % rotate particle motions into body-frame
    v_b = jac' * [ waves.vx; waves.vz; waves.v_theta ];
    a_b = jac' * [ waves.ax; waves.az; waves.a_theta ];

    % determine control actions
    [ tau, mu ] = cascaded_pid_control( brov2, jac, ...
        brov2.ekf.p_error, brov2.ekf.i_error, brov2.ekf.d_error, dt );
        
    brov2.plot.mu(k,:) = mu(:,1);

    % estimate state with an extended kalman filter
    brov2 = state_estimator( brov2, tau(1:2), v_b(1:2), a_b(1:2), dt, k );

    % simulate rov dynamics, calculate state errors
    brov2 = dynamics( brov2, v_b, a_b, tau, jac, dt );
    brov2.plot.state(k,:) = brov2.state';
    [ brov2 ] = calc_errors_ekf( brov2, brov2.p_error, brov2.ekf.p_error, k );
    
    fprintf('timstep: %i \n', k);

end

%% produce gif to animate the rov behaviour under waves

figure('Position', [25 200 600 600]);
filename = 'test.gif';

x = 48:0.1:52;

for j = 1:length(brov2.plot.state(:,1))

    eta = zeros(numel(x),1);

    for jj = 1:numel(eta)
        for ii = 1:numel(waves.A)
            eta(jj) = eta(jj) + waves.A(ii) * cos( waves.k(ii)*x(jj) - ...
                waves.w(ii)*t(j) + waves.E(ii) );
        end
    end

    area(x, eta, -15, 'facecolor', [0.3010 0.7450 0.9330]);
    ylim([-10 3]);
    xlim([min(x) max(x)]);
    hold on

    h = rectangle('Position', [brov2.plot.state(j,1)-brov2.length/2 brov2.plot.state(j,2)-brov2.height*2 ...
            brov2.length brov2.height*4], 'Curvature', 0.25);
    h.FaceColor = 'r';
% %     g = hgtransform('Matrix', makehgtform('yrotate', brov2.plot.state(j,3)));
% %     h.Parent = g;

    scatter(brov2.plot.state(j,1), brov2.plot.state(j,2), 50, 'ok', 'filled')
    plot(ref_pose(1), ref_pose(2), 'x', 'Linewidth', 2);
    yline(0, 'r', 'linestyle', '--', 'linewidth', 1.5);
    hold off

    xlabel('x (m)', 'Interpreter', 'latex');
    ylabel('z (m)', 'Interpreter', 'latex');
    legend('Wave', 'ROV', 'Target', 'Still Water Line', 'Interpreter', 'latex');
    title("t = " + t(j) + "s", 'Interpreter', 'latex');

    drawnow

    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

end

%% plot simulation outputs

figure;
subplot(3,1,1);
plot(t, brov2.plot.state(:,1));
hold on;
% plot(t, brov2.ekf.plot.y(:,1), '.');
% plot(t, brov2.ekf.plot.x_hat(:,1));

hold off;

subplot(3,1,2);
plot(t, brov2.plot.state(:,2));
hold on;
% plot(t, brov2.ekf.plot.y(:,2), '.');
% plot(t, brov2.ekf.plot.x_hat(:,2));
hold off;

subplot(3,1,3);
plot(t, rad2deg(brov2.plot.state(:,3)));


