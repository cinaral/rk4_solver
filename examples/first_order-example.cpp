#include "rk4_solver/loop.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e3;
constexpr Real_T h = 1. / sample_freq;
constexpr Real_T t0 = 0.;
constexpr Real_T tf = 10.;
constexpr size_t t_dim = sample_freq * (tf - t0) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x0[x_dim] = {1.};
constexpr Real_T a_constant = 1;

struct Dynamics {
	/*
	 * first order ode:
	 * dt_x = a*x
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = a_constant * x[0];
	}
};
Dynamics dynamics;

int
main()
{
	Real_T t;
	Real_T x[x_dim];
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, h);
	//* integration loop
	rk4_solver::loop<t_dim>(integrator, t0, x0, &t, x);

	return 0;
}
