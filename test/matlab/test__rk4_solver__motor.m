%********
%* setup
%********
io_config;
test_name = 'test-rk4_solver-motor';
prefix = append(dat_dir, '/', test_name, '-');
exe_name = append(test_name, '.exe');
t_arr_fname = 't_arr.dat';
x_arr_fname = 'x_arr.dat';
u_arr_fname = 'u_arr.dat';
x_arr_chk_fname = 'x_arr_chk.dat';

t_dim = 1e3;
x_dim = 3;
error_thres = 1e-9;
f = 5;

%*******************
%* create test data 
%*******************
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

t_arr = linspace(0, 1, t_dim).';
u_arr = sin(t_arr*2*pi*f);

%! need x_arr_chk
%sys = idss(A, B, C, D);
%y = sim(sys, u_arr);

%************************************
%* write input (for test executable)
%************************************
writematrix(t_arr, append(prefix, t_arr_fname), 'Delimiter', delimiter);  
writematrix(u_arr, append(prefix, u_arr_fname), 'Delimiter', delimiter);  
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
x_arr = readmatrix(append(prefix, x_arr_fname));

%*********
%* verify
%*********
max_error = max(vecnorm(x_arr - x_arr_chk, 2, 2));
mean_error = mean(vecnorm(x_arr - x_arr_chk, 2, 2));

if max_error < error_thres
    disp(append(test_name, '	ok'));
else
    disp(append(test_name, '	fail'));
end

figure('Name', 'x');
hold on;
plot(t_arr, x_arr);
plot(t_arr, x_arr_chk, '--');