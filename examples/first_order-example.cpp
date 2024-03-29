#include "rk4_solver/loop.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 10.;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x_init[x_dim] = {1.};
constexpr Real_T a_constant = 1;

struct Dynamics {
	/*
	 * first order ode:
	 * dt_x = a*x
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = a_constant * x[0];
	}
};
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);

	Real_T t;
	Real_T x[x_dim];
	rk4_solver::loop<t_dim>(integrator, t_init, x_init, t, x);

	return 0;
}
