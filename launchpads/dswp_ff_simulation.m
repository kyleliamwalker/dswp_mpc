% Author: Dr. Kyle L. Walker
% Description: Launchpad for dynamically simulating a PID controlled ROV in
% 3DOF (surge, heave, pitch) under the influence of wave disturbances, 
% incorporating feedforward disturbance rejection. For the control, 
% disturbances are estimated using a Determistic Sea Wave Predictor (DSWP) 
% and a Linear Wave Theory (LWT) based disturbance model. Desired set-point 
% is altered by modifying ref_pose. An Extended Kalman Filter is employed 
% for estimating the vehicle state. A GIF is produced at the end of the
% simulation to visualise the vehicle motion.

%%

clear variables
%close all
clc

% set up time parameters
dt = 0.1;
t = dt:dt:500;
% set prediction distance and measurement time for dswp
dx = 50;
t_measure = 300;

% set initial and desired rov pose
robot_pose = [ dx; -5; 0 ];
ref_pose = [ dx; -5; 0 ];

%% load wave parameters and vehicle specification
waves = dswp_load_waves( t, robot_pose(1), dx );
% waves = dswp_load_2nd_order_waves( t, robot_pose(1), dx );

brov2 = load_robot( robot_pose, ref_pose, t );
% initialise array to store wave height measurements
m_wave = zeros(1, numel(t));

%% run dynamic simulation
for k = 1:numel(t)
    
    % store measure wave height at each timestep and add white noise
    m_wave(k) = awgn(waves.eta(k), 10);
    
    if t(k) > t_measure

        % jacobian for earth/body frame rotation
        jac = [ cos(brov2.state(3))  sin(brov2.state(3))  0;
                -sin(brov2.state(3)) cos(brov2.state(3))  0;
                    0           0       1];

        %filter measurements to remove noise
        filt_wave = smoothdata(m_wave(1:k), 'gaussian', 7);
        % obtain wave spectrum with dswp
        [ t_pred, dswp_waves ] = dswp(t(k), dx, filt_wave ) ;
        dswp_waves.d = waves.d;
        dswp_waves.rho = waves.rho;

        % calculate estimated particle motions 
        [ dswp_waves ] = evaluate_particles( t_pred, brov2.ekf.x_hat(1), ...
              brov2.ekf.x_hat(2), dswp_waves );
%         [ dswp_waves ] = evaluate_2nd_order_particles( t_pred, brov2.ekf.x_hat(1), ...
%             brov2.ekf.x_hat(2), dswp_waves );

        % estimate wave profile along prediction horizon
        pred_wave = zeros(1,numel(t_pred));
        for i = 1:numel(dswp_waves.T)   
            pred_wave = pred_wave + dswp_waves.H(i) / 2 * ...
                cos(dswp_waves.k(i)*dx - dswp_waves.w(i)*t_pred + ...
                dswp_waves.E(i));
        end

        % use particles from dswp to get disturbances
        dist = ff_disturbances( brov2, dswp_waves, jac, t, k-1 );

        % use dswp prediction to calculate ff control and allocate
        [ tau, mu ] = ff_control( brov2, jac, brov2.ekf.p_error, ...
            brov2.ekf.i_error, brov2.ekf.d_error, dist, dt );
        brov2.plot.mu(k,:) = mu(:,1);
        brov2.plot.tau(k,:) = tau(:,1);
        
        % get actual wave particles for simulation
        [ waves ] = evaluate_particles( t, brov2.state(1), ...
            brov2.state(2), waves );
%         [ waves ] = evaluate_2nd_order_particles( t, brov2.state(1), ...
%             brov2.state(2), waves );
        [ waves ] = pitch_particles( brov2, waves, t(k), jac );
        waves.plot.vx(k) = waves.vx(k);    waves.plot.ax(k) = waves.ax(k);
        waves.plot.vz(k) = waves.vz(k);    waves.plot.az(k) = waves.az(k);
        waves.plot.v_theta(k) = waves.v_theta;   waves.plot.a_theta(k) = waves.a_theta;

        % rotate into body frame
        v_b = jac' * [ waves.vx(k); waves.vz(k); waves.v_theta ];
        a_b = jac' * [ waves.ax(k); waves.az(k); waves.a_theta ];
        
        % estimate state
        brov2 = state_estimator( brov2, tau(1:2), v_b(1:2), a_b(1:2), dt, k );
        % simulate rov dynamics
        brov2 = dynamics( brov2, v_b, a_b, tau, jac, dt );
        
    end
    
    % plot everything and calculate errors
    brov2.plot.state(k,:) = brov2.state';
    [ brov2 ] = calc_errors_ekf( brov2, brov2.p_error, brov2.ekf.p_error, k );
    
    fprintf('timstep: %i \n', k);
end

%% produce gif to animate the rov behaviour under waves

figure('Position', [25 200 600 600]);
filename = 'test.gif';

x = 48:0.1:52;

for j = t_measure/dt:length(brov2.plot.state(:,1))

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
    if j == t_measure/dt
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

end


%% plot simulation outputs

figure;
subplot(4,1,1);
plot(t, waves.pred);
xlim([t_measure k*dt]);
subplot(4,1,2);
plot(t, brov2.plot.state(:,1)-50, 'r');
xlim([t_measure k*dt]);
subplot(4,1,3);
plot(t, brov2.plot.state(:,2), 'r');
xlim([t_measure k*dt]);
subplot(4,1,4);
plot(t, rad2deg(brov2.plot.state(:,3)), 'r');
xlim([t_measure k*dt]);

