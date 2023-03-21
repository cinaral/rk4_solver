%* setup
%*
%* See https://www.mathworks.com/help/matlab/math/ode-event-location.html
%* for the official MATLAB example counterpart of what we are trying to achieve in this test.
%* 
%* This is very similar ballode but it is fixed-step integration here.
%*
bin_dir = '../../build/bin';
dat_dir = '../../build/dat';
ref_dat_dir = "../reference_dat";
test_name = 'ball-test';
dat_prefix = append(dat_dir, '/', test_name, '-');
ref_dat_prefix = append(ref_dat_dir, '/', test_name, '-');
exe_name = append(test_name, '.exe');
t_arr_fname = 't_arr.dat';
x_arr_fname = 'x_arr.dat';
x_arr_ref_fname = 'x_arr_ref.dat';

is_drawing = true;
is_single_precision = false;
sample_freq = 1e4;
t_init = 0;
t_final = 2;
t_dim = sample_freq*(t_final - t_init) + 1;
x_dim = 2;
u_dim = 1;
gravity_const = 9.806;
e_restitution = 0.75;
x_init = [1; 0];

if is_single_precision
	error_thres = 1e-3;
else 
	error_thres = 1e-3;
end

%* create reference data 
ball_system = @(~, x) [x(2); -gravity_const];

refine = 16;
options = odeset('Events', @event, 'RelTol', 1e-10, 'AbsTol', 1e-12); %* 'RelTol', 1e-10, 'AbsTol', 1e-12

t_arr_ref = linspace(t_init, t_final, t_dim).';
x_arr_ref = zeros(t_dim, x_dim);

idx_begin = 1;
t_arr_part = t_arr_ref;
t = 0;

while t < t_final && ~isempty(t_arr_part)
	%* Numerical integration
	[t_arr_part, x_arr_part] = ode45(ball_system, t_arr_part, x_init, options);
	t = t_arr_part(end);
	len = length(t_arr_part);
	x_init(1) = 0;
	x_init(2) = -e_restitution * x_arr_part(end, 2);

	%* Accumulate output
	x_arr_ref(idx_begin:idx_begin + len - 1, :) = x_arr_part;
	idx_begin = idx_begin + len;
	t_arr_part = t_arr_ref(idx_begin:end);
end

writematrix(x_arr_ref, append(ref_dat_prefix, x_arr_ref_fname));  

disp(append('Created reference data for ', test_name));

prev_pwd = pwd;
cd(bin_dir);
if isfile(exe_name)
	%* call the test executable
	if system(exe_name) > 0
		warning(append(bin_dir, '/', exe_name, ' has returned failure.'));
	end
	
	%* read the results
	t_arr = readmatrix(append(dat_prefix, t_arr_fname));
	x_arr = readmatrix(append(dat_prefix, x_arr_fname));

	%* verify
	max_error = max(vecnorm(x_arr(:, 1) - x_arr_ref(:, 1), 1, 2));
	mean_error = mean(vecnorm(x_arr(:, 1) - x_arr_ref(:, 1), 1, 2));

	if max_error < error_thres
		disp(append(test_name, '	ok'));
	else
		disp(append(test_name, '	fail'));
	end

	if is_drawing
		figure('Name', 'x');
		hold on;
		plot(t_arr, x_arr(:, 1));
		plot(t_arr_ref, x_arr_ref(:, 1), '--');
	end
else
	error(append(bin_dir, '/', exe_name, ' does not exist. Use CMake to build the test.'));
end
cd(prev_pwd);


function [value, isterminal, direction] = event(~, x)
	value = x(1); % Detect height = 0
	isterminal = 1; % Stop the integration
	direction = -1; % Negative direction only
end
