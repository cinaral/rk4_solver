%* setup
addpath('../');
test_name = 'test-rk4_solver-exact_integral';
is_drawing = false;
error_thres = 1e-10;
%* Exact Integral
%* dt__y = t.^3 - t.^2

sample_freq = 1e3;
time_step = 1/sample_freq;
t_init = 0;
t_final = 20;
t_dim = sample_freq*(t_final - t_init) + 1;
x_dim = 1;
y_fun = @(t) t.^4/4 - t.^3/3;
dt__y_fun = @(t) t.^3 - t.^2;

t_arr = linspace(t_init, t_final, t_dim).';
y_arr = zeros(t_dim, x_dim);

%* call
for i = 1:t_dim - 1
	t = t_arr(i, :).';
	y = y_arr(i, :).';
	h = t_arr(i + 1) - t;
	f = @(t, x) dt__y_fun(t);
	y_arr(i + 1, :) = step_rk4(t, y, h, f).';
end

%* verify
max_error = max(abs(y_arr - y_fun(t_arr)));
mean_error = mean(abs(y_arr - y_fun(t_arr)));

if max_error < error_thres
	disp(append(test_name, '	ok'))
else
	disp(append(test_name, '	fail'))
end

if is_drawing
    figure('Name', 'y')
    plot(t_arr, y_arr)
    hold on
    plot(t_arr, y_fun(t_arr), '--')
    xlabel('t (s)');
    ylabel('$\mathbf{y}(t)$');
    legend('$\mathbf{y}(t)$ integrated by step\_rk4', '$\mathbf{y}(t)$', 'Location', 'best')

    figure('Name', 'y error');
    hold on;
    plot(t_arr, y_arr - y_fun(t_arr));
end
rmpath('../');
