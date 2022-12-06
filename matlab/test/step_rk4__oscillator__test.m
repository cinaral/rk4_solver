%* setup
addpath('../');
test_name = 'rk4_solver-oscillator-test';
is_drawing = false;
error_thres = 2e-3;
%* Linear 1-DoF oscillator excited by sine
%* dtdt__y + 2*zeta*w_n*dt__y + w_n^2*y = sin(t)
%* x = [y; dt__y]


sample_freq = 1e4;
time_step = 1/sample_freq;
t_init = 0;
t_final = 20;
t_dim = sample_freq*(t_final - t_init) + 1;
x_dim = 2;
u_dim = 1;
zeta = 0.1;
w_n = 1;
A = [0, 1; -w_n^2, -2*zeta*w_n];
B = [0; 1];
dt__x_fun = @(x, u) A*x + B*u;
x_init = [1; 0];
w_u = 1;

t_arr = linspace(t_init, t_final, t_dim).';
u_arr = sin(pi*w_u*t_arr);
x_arr = [x_init.'; zeros(t_dim - 1, x_dim)];

%* call
for i = 1:t_dim - 1
	t = t_arr(i, :).';
	x = x_arr(i, :).';
	h = t_arr(i + 1) - t;
	u_0 = u_arr(i, :).';
	f_0 = @(t, x) dt__x_fun(x, u_0);
	x_arr(i + 1, :) = step_rk4(t, x, h, f_0).';
end

%* verify
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
