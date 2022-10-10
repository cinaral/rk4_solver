%* setup
addpath('../');
test_name = 'test-rk4_solver-oscillator';
is_drawing = false;
error_thres = 2e-3;
%* Linear 1-DoF oscillator excited by sine
%* dtdt__y + 2*zeta*w_n*dt__y + w_n^2*y = sin(t)
%* x = [y; dt__y]

time_step = 1e-4;
t_init    = 0;
t_final   = 20;
t_arr     = (t_init:time_step:t_final).';
t_arr_len = size(t_arr, 1);
x_dim     = 2;
u_dim     = 1;
zeta      = 0.1;
w_n       = 1;
A         = [0, 1; -w_n^2, -2*zeta*w_n];
B         = [0; 1];
dt__x_fun = @(x, u) A*x + B*u;
x_init    = [1; 0];
w_u       = 1;
u_arr     = sin(pi*w_u*t_arr);
x_arr     = [x_init.'; zeros(t_arr_len - 1, x_dim)];

%* call
for i = 1:t_arr_len - 1
	t = t_arr(i, :).';
	x = x_arr(i, :).';
	h = t_arr(i + 1) - t;
	u_0 = u_arr(i, :).';
	f_0 = @(t, x) dt__x_fun(x, u_0);
	x_arr(i + 1, :) = step_rk4(t, x, h, f_0).';
end

% verify
ode_fun  = @(t, x) A*x + B*sin(pi*w_u*t); 
[ode.t_arr, ode.x_arr] = ode45(ode_fun, t_arr, x_init); 

max_error = max(vecnorm(x_arr - ode.x_arr, 2, 2));
mean_error = mean(vecnorm(x_arr - ode.x_arr, 2, 2));

if max_error < error_thres
	disp(append(test_name, '	ok'));
else
	disp(append(test_name, '	fail'));
end

if is_drawing
    figure('Name', 'x')
    plot(t_arr, x_arr)
    hold on
    plot(ode.t_arr, ode.x_arr, '--')
    xlabel('t (s)');
    ylabel('$\mathbf{x}(t)$');
    legend('$\mathbf{x}_1(t)$ integrated by step\_rk4', '$\mathbf{x}_2(t)$ integrated by step\_rk4', ...
	    '$\mathbf{x}_1(t)$ integrated by ode45', '$\mathbf{x}_2(t)$ integrated by ode45', 'Location', 'best')
end
rmpath('../');
