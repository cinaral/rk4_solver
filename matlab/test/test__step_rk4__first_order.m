%* setup
addpath('../');
test_name = 'test-rk4_solver-first_order';
is_drawing = false;
error_thres = 1e-12;
%* Exact Integral
%* f(t, x) = a*x
%* x = exp(a*t)

time_step = 1e-3;
t_init = 0;
t_final = 1;
t_arr = (t_init:time_step:t_final).';
t_arr_len = size(t_arr, 1);
a = 1;
x0 = 1;
x_fun = @(t) exp(a*t);
dt__x_fun = @(t, x) a*x;
x_arr = zeros(t_arr_len, 1);
x_arr(1, :) = 1;

%* call
for i = 1:t_arr_len - 1
	t = t_arr(i, :).';
	y = x_arr(i, :).';
	h = t_arr(i + 1) - t;
	f = @(t, x) dt__x_fun(t, x);
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
