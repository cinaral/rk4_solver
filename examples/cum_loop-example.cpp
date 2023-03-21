#include "rk4_solver/cum_loop.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e3;
constexpr Real_T h = 1. / sample_freq;
constexpr Real_T t0 = 0.;
constexpr Real_T tf = 1.;
constexpr size_t t_dim = sample_freq * (tf - t0) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x0[x_dim] = {1.};

struct Dynamics {
	/*
	 * dt_x = f(t, x) = t, x(0) = 0
	 * x = 1/2*t^2
	 */
	void
	ode_fun(const Real_T t, const Real_T (&)[x_dim], const size_t, Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = t;
	}
};
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, h);
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	//* integration loop with cumulatively saved data arrays
	rk4_solver::cum_loop<t_dim>(integrator, t0, x0, t_arr, x_arr);

	return 0;
}
