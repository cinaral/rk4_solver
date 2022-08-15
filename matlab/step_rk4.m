function [x_new, t_new] = step_rk4(t, x, h, f_0, f_1, f_2)
%* Runge-kutta 4th order integration step for functions with time-dependent parameters.
%* step_rk4(t, x, h, f_0, f_1, f_2)
%* t   - time
%* x   - state 
%* h   - time step
%* f_0 - integrand at t
%* f_1 - integrand at t + h/2
%* f_2 - integrand at t + h
%*
%* step_rk4(t, x, h, f_0): If there are no time-dependent parameters then 
%* f_1 and f_2 can be skipped.
%* 
%* step_rk4(t, x, h, f_0, f_1, f_2): When time-dependent parameters are present, 
%* it is preferred that parametrized functions f_1 and f_2 are used instead of
%* embedding interpolating functions into f_0.
%*
%* dt__x     = f(t, x, a(t))
%* f_0(t, x) = f(t, x, a(t))
%* f_1(t, x) = f(t, x, a(t + h/2))
%* f_2(t, x) = f(t, x, a(t + h))
%*
%* Check step_rk4_test.m for usage and verification.
%*
%* cinaral 2022-01-21

if nargin < 5
	f_1 = f_0;
	f_2 = f_0;
end

if nargin == 5
	if coder.target('MATLAB')
		st = dbstack;
		error('%s: You must specify both f_1 and f_2 to use parametrized integrands.\n', st(1).name);
	else
		error('You must specify both f_1 and f_2 to use parametrized integrands.\n')
	end
end

k1 = f_0(t,       x);
k2 = f_1(t + h/2, x + h*k1/2);
k3 = f_1(t + h/2, x + h*k2/2);
k4 = f_2(t + h,   x + h*k3);

x_new = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
t_new = t + h;