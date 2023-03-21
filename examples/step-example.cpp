#include "rk4_solver/integrator.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t x_dim = 2;
constexpr Real_T time_step = 1e-3;
constexpr Real_T t_init = 0;
constexpr Real_T x_init[x_dim] = {0, 0};
constexpr Real_T u = 1e3;

struct Dynamics {
	/*
	 * dt_x = f(t, x) = [x[1]; u], x(0) = [0; 0];
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = x[1];
		dt_x[1] = u;
	}
};
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	Real_T t;
	Real_T x[x_dim];
	integrator.step(t_init, x_init, t, x); //* integration step

	return 0;
}
