%********
%* setup
%********
io_config;
test_name = 'test-rk4_solver-sine';
prefix = append(dat_dir, '/', test_name, '-');
exe_name = append(test_name, '.exe');
t_arr_fname = 't_arr.dat';
x_arr_fname = 'x_arr.dat';
x_arr_chk_fname = 'x_arr_chk.dat';

t_dim = 1e3;
x_dim = 1;
error_thres = 1e-9;
f = 5;


%writematrix(t_arr, append(prefix, t_arr_fname), 'Delimiter', delimiter);  
%writematrix(x_arr_chk, append(prefix, x_arr_chk_fname), 'Delimiter', delimiter);  

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

%*******************
%* create test data 
%*******************
x_arr_chk = sin(t_arr*2*pi*f);

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

%figure('Name', 'x');
%hold on;
%plot(t_arr, x_arr);
%plot(t_arr, x_arr_chk, '--');