%********
%* setup
%********
addpath('..\..\matlab');

io_config;
test_name = 'test-rk4_solver-ball';
prefix = append(dat_dir, '/', test_name, '-');
exe_name = append(test_name, '.exe');
t_arr_fname = 't_arr.dat';
x_arr_fname = 'x_arr.dat';
x_arr_chk_fname = 'x_arr_chk.dat';

t_dim = 1e5;
x_dim = 2;
u_dim = 1;
error_thres = 1e-4;
%* __USE_SINGLE_PRECISION__
%error_thres = 1e-3;
g = 9.806;
e = 0.75;
x0 = [1; 0];
tf = 2;

%*******************
%* create test data 
%*******************

ball_system = @(~, x) [x(2); -g];

refine = 16;
options = odeset('Events', @event, 'Refine', refine);

t_arr_chk = linspace(0, tf, t_dim).';
x_arr_chk = zeros(t_dim, x_dim);

idx_begin = 1;
t_arr_part = t_arr_chk;
t = 0;

while t < tf && ~isempty(t_arr_part)
	% Numerical integration
	[t_arr_part, x_arr_part] = ode45(ball_system, t_arr_part, x0, options);

% 	i = length(t_arr_part);
	t = t_arr_part(end);
	len = length(t_arr_part);
	x0(1) = 0;
	x0(2) = -e * x_arr_part(end, 2);
% 	x_arr_part(end, :) = x0;

	% Accumulate output
	x_arr_chk(idx_begin:idx_begin + len - 1, :) = x_arr_part;
	
	idx_begin = idx_begin + len;
	t_arr_part = t_arr_chk(idx_begin:end);
end

%figure('Name', 'u');
%plot(t_arr, u_arr)
%figure('Name', 'x_arr_chk');
%plot(t_arr, x_arr_chk(:, 1))

%************************************
%* write input (for test executable)
%************************************
writematrix(x_arr_chk, append(prefix, x_arr_chk_fname), 'Delimiter', delimiter);  

%***************************
%* call the test executable
%***************************
prev_pwd = pwd;
cd(bin_dir);

if ~isfile(exe_name)
	error(append(bin_dir, '/', exe_name, ' does not exist. Use CMake to build the test.'));
end

if system(exe_name) > 0
	warning(append(bin_dir, '/', exe_name, ' has returned failure.'));
end

cd(prev_pwd);

%******************************************
%* read output (created by the executable)
%******************************************
t_arr = readmatrix(append(prefix, t_arr_fname));
x_arr = readmatrix(append(prefix, x_arr_fname));

%*********
%* verify
%*********
max_error = max(vecnorm(x_arr(:, 1) - x_arr_chk(:, 1), 1, 2));
mean_error = mean(vecnorm(x_arr(:, 1) - x_arr_chk(:, 1), 1, 2));

if max_error < error_thres
    disp(append(test_name, '	ok'));
else
    disp(append(test_name, '	fail'));
end

%figure('Name', 'x');
%hold on;
%plot(t_arr, x_arr(:, 1));
%plot(t_arr_chk, x_arr_chk(:, 1), '--');


function [value, isterminal, direction] = event(~, x)
	value = x(1); % Detect height = 0
	isterminal = 1; % Stop the integration
	direction = -1; % Negative direction only
end
