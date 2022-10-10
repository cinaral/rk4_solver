function [x_new, t_new] = step_rk4(t, x, h, f)
%* Integrates dt__x = ode_fun(t, x(t)), x(0) = x_0 for one time-step using Runge-Kutta 4th Order Method.
%* 
%* [OUT:x_new, OUT:t_new] = step_rk4(t, x, h, f):
%* IN:
%* t - time [s]
%* x - [X_DIM] state 
%* h - time step [s]
%* f - integrand at t
%* OUT:
%* x_new - [X_DIM] new state
%* t_new - new time [s]
%*
%* cinaral 2022-10-10

k1 = f(t, x);
k2 = f(t + h/2, x + h*k1/2);
k3 = f(t + h/2, x + h*k2/2);
k4 = f(t + h, x + h*k3);

x_new = x + h/6*(k1 + 2*k2 + 2*k3 + k4);
t_new = t + h;