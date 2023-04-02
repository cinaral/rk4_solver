#include "rk4_solver/loop.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 1.;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x_init[x_dim] = {1.};

struct Dynamics {
	/*
	 * dt_x = f(t, x) = t, x(0) = 0
	 * x = 1/2*t^2
	 */
	void
	ode_fun(const Real_T t, const Real_T (&)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = t;
	}
};
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim][x_dim];
	//* save the intermediate values
	rk4_solver::loop(integrator, t_init, x_init, t_arr, x_arr);

	return 0;
}
