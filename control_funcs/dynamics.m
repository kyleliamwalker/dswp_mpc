function [ robot ] = dynamics( robot, v_b, a_b, tau, jac, dt )

u = robot.state(4); w = robot.state(5); q = robot.state(6);
theta = robot.state(3); 

M_RB = robot.M_RB;
% leave this in?
M_A = robot.M_A;
D_v_quad = robot.D_v_quad;
D_v_lin = robot.D_v_lin;
r_g = robot.r_g;
r_b = robot.r_b;
W = robot.W;
B = robot.B;

% coriolis tests
M = blkdiag(robot.mDry,0, robot.mDry, 0, robot.Iyy, 0);
M2 = [ robot.Xu_dot, 0, robot.Xw_dot, 0, robot.Xq_dot, 0;
       0, 0, 0, 0, 0, 0;
       robot.Zu_dot, 0, robot.Zw_dot, 0, robot.Zq_dot, 0;
       0, 0, 0, 0, 0, 0;
       robot.Mu_dot, 0, robot.Mw_dot, 0, robot.Mq_dot, 0;
       0, 0, 0, 0, 0, 0];
   
C_RB = m2c(M, [ u 0 w 0 q 0 ]'); 
C_RB([2 4 6],:) = [];
C_RB(:,[2 4 6]) = [];

C_A  = m2c(M2, [ u 0 w 0 q 0 ]');
C_A([2 4 6],:) = [];
C_A(:,[2 4 6]) = [];

% assume weight and buoyancy are equal
g_nu = [ (0) * sin(theta)
        -(0) * cos(theta) * cos(0)
        (r_g(2)*W-r_b(2)*B) * sin(theta) ];

g_lin = [ (0) * sin(theta)
        -(0) * cos(theta) * cos(0)
        (r_g(2)*W-r_b(2)*B) * theta ];
    
robot.dist = M_A*a_b + ( D_v_lin + D_v_quad .* abs([u;w;q]-v_b) ) * ([u;w;q]-v_b);

nu_dot = @(t,nu_in) ...
    -(M_RB + M_A) \ ...
    (M_A*a_b + (C_RB + C_A)*(nu_in) + ( D_v_lin + D_v_quad .* abs(nu_in-v_b) ) * (nu_in-v_b) ...
    + g_nu -  tau) ;

%linear dynamics

% nu_dot = @(t,nu_in) ...
%     -(M_RB + M_A) \ (M_A*a_b + ( D_v_lin ) * (nu_in-v_b) ...
%     + g_lin -  tau) ;
% + D_v_quad .* abs(nu_in-v_b) 
    
[ t, nu ] = ode45( nu_dot, [0 dt], [ u w q ]' );

nudot = nu_dot( t(end), nu(end,:)' );

udot = nudot(1);    wdot = nudot(2);    qdot = nudot(3);
u = nu(end, 1);     w = nu(end, 2);     q = nu(end, 3);

% transform into world frame
eta_dot = jac * nu';

x = odeDisplacement( robot.state(1), eta_dot(1,:)', t );
z = odeDisplacement( robot.state(2), eta_dot(2,:)', t );
theta = odeDisplacement( robot.state(3), eta_dot(3,:)', t );

robot.state = [ x; z; theta; u; w; q; udot; wdot; qdot ];

end

