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
	//* first order ode:
	//* dt__x = a*x
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = a_constant * x[0];
	}
};
Dynamics dyn;

int
main()
{
	Real_T t;
	Real_T x[x_dim];
	//* integration loop
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t0, x0, h, &t, x);

	return 0;
}
