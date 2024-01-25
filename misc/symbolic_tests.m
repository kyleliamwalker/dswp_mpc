
clear all

syms x xdot xddot x0 M M_A D_l D_q G tau a_b v_b nu

f = (M+M_A)^-1 * ( M_A*a_b + ( D_l + D_q.*abs( nu ) ) * ( nu )) ;

f_lin = subs(f, x0) + (nu-x0)*subs(diff(f,nu),x0);
f_lin = simplify(f_lin);

disp(f_lin)
% ezplot(f,[0,5])
% hold on
% ezplot(subs(f_lin,x0,2),[1 3])
% legend('f(x)', 'f_lin(x)')

%% multi-directional particles
% need to adapt this to include effect of varying sea bed depth

clear all

syms A w k z d x y z t E dir

X = x*cos(dir) + y*sin(dir);
%X = x;

pot = A*w/k * (cosh(k*(z+d))/cosh(k*d)) * sin(k*X - w*t + E);

vx = simplify(diff(pot, x));
vy = simplify(diff(pot, y));
vz = simplify(diff(pot, z));

ax = simplify(diff(vx, t));
ay = simplify(diff(vy, t));
az = simplify(diff(vz, t));

%% create functions

matlabFunction(vx, 'vars', {A, w, k, E, dir, d, x, y, z, t}, 'file', 'wave_funcs/md_vx');
matlabFunction(vy, 'vars', {A, w, k, E, dir, d, x, y, z, t}, 'file', 'wave_funcs/md_vy');
matlabFunction(vz, 'vars', {A, w, k, E, dir, d, x, y, z, t}, 'file', 'wave_funcs/md_vz');
matlabFunction(ax, 'vars', {A, w, k, E, dir, d, x, y, z, t}, 'file', 'wave_funcs/md_ax');
matlabFunction(ay, 'vars', {A, w, k, E, dir, d, x, y, z, t}, 'file', 'wave_funcs/md_ay');
matlabFunction(az, 'vars', {A, w, k, E, dir, d, x, y, z, t}, 'file', 'wave_funcs/md_az');