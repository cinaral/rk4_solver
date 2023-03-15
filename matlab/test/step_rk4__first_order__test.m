%* setup
addpath('../');
test_name = 'rk4_solver-first_order-test';
is_drawing = false;
error_thres = 1e-12;
%* Exact Integral
%* f(t, x) = a*x
%* x = exp(a*t)


sample_freq = 1e3;
time_step = 1/sample_freq;
t_init = 0;
t_final = 1;
t_dim = sample_freq*(t_final - t_init) + 1;
x_dim = 1;
a_constant = 1;
x0 = 1;
x_fun = @(t) exp(a_constant*t);
dt_x_fun = @(t, x) a_constant*x;

t_arr = linspace(t_init, t_final, t_dim).';
x_arr = zeros(t_dim, x_dim);
x_arr(1, :) = x0;

%* call
for i = 1:t_dim - 1
	t = t_arr(i, :).';
	y = x_arr(i, :).';
	h = t_arr(i + 1) - t;
	f = @(t, x) dt_x_fun(t, x);
	x_arr(i + 1, :) = step_rk4(t, y, h, f).';
end

%* verify
max_error = max(abs(x_arr - x_fun(t_arr)));
mean_error = mean(abs(x_arr - x_fun(t_arr)));

if max_error < error_thres
	disp(append(test_name, '	ok'))
else
	disp(append(test_name, '	fail'))
end

if is_drawing
    figure('Name', 'x')
    plot(t_arr, x_arr)
    hold on
    plot(t_arr, x_fun(t_arr), '--')
    xlabel('t (s)');
    ylabel('$\mathbf{x}(t)$');
    legend('$\mathbf{x}(t)$ integrated by step\_rk4', '$\mathbf{x}(t)$', 'Location', 'best')

    figure('Name', 'x error');
    hold on;
    plot(t_arr, x_arr - x_fun(t_arr));
end
rmpath('../');
