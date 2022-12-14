%* setup
addpath('../../matlab');
bin_dir = '../../build/bin';
dat_dir = '../../build/dat';
ref_dat_dir = "../reference_dat";
test_name = 'motor-test';
dat_prefix = append(dat_dir, '/', test_name, '-');
ref_dat_prefix = append(ref_dat_dir, '/', test_name, '-');
exe_name = append(test_name, '.exe');
t_arr_fname = 't_arr.dat';
x_arr_fname = 'x_arr.dat';
u_arr_fname = 'u_arr.dat';
x_arr_chk_fname = 'x_arr_chk.dat';

is_drawing = false;
is_single_precision = true;
sample_freq = 1e3;
t_init = 0;
t_final = 1;
t_dim = sample_freq*(t_final - t_init) + 1;
x_dim = 3;
u_dim = 1;

if is_single_precision
	error_thres = 1e-3;
else 
	error_thres = 1e-4;
end

%* create reference data 
f = 1000; %* input frequency
R = 1.4; %* [ohm]
L = 1.7e-3; %* [ohm s] 
J = 1.29e-4; %* [kg m-2]
b = 3.92e-4; %* [N m s]
K_t = 6.4e-2; %* [N m A-1]
K_b = 6.4e-2; %* [V s]

A = [0, 1, 0; 
	 0, -b/J,   K_t/J; 
	 0, -K_b/L, -R/L];
B = [0; 
	 0; 
	 1/L];
C = [1, 0, 0];
D = 0;

t_arr_chk = linspace(t_init, t_final, t_dim).';
motor_system = @(x, u) A*x + B*u;

x_arr_chk = zeros(t_dim, x_dim);
u_arr = zeros(t_dim, u_dim);

for i = 1:t_dim 
	t = t_arr_chk(i, :).';
	x = x_arr_chk(i, :).';
	u = sin(t_arr_chk(i)*2*pi*f);

	if i < t_dim
		h = t_arr_chk(i + 1) - t;	
		x_arr_chk(i + 1, :) = step_rk4(t, x, h, @(t, x) motor_system(x, u)).';
	end
	u_arr(i, :) = u.';
end

writematrix(u_arr, append(ref_dat_prefix, u_arr_fname));  
writematrix(x_arr_chk, append(ref_dat_prefix, x_arr_chk_fname));  

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
	max_error = max(vecnorm(x_arr - x_arr_chk, 2, 2));
	mean_error = mean(vecnorm(x_arr - x_arr_chk, 2, 2));

	if max_error < error_thres
		disp(append(test_name, '	ok'));
	else
		disp(append(test_name, '	fail'));
	end

	if is_drawing
		figure('Name', 'x');
		hold on;
		plot(t_arr, x_arr(:, 1));
		plot(t_arr_chk, x_arr_chk(:, 1), '--');
	end
else
	error(append(bin_dir, '/', exe_name, ' does not exist. Use CMake to build the test.'));
end
cd(prev_pwd);
rmpath('../../matlab');