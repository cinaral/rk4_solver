#include "rk4_solver/loop.hpp"

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;

constexpr uint_t sample_freq = 1e3;
constexpr real_t time_step = 1. / sample_freq;
constexpr real_t t_init = 0.;
constexpr real_t t_final = 10.;
constexpr uint_t t_dim = sample_freq*(t_final - t_init) + 1;
constexpr uint_t x_dim = 1;
constexpr real_t x_init[x_dim] = {1.};
constexpr real_t a_constant = 1;

struct Dynamics {
	//* first order ode:
	//* dt__x = a*x
	void
	ode_fun(const real_t, const real_t (&x)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
	{
		dt__x[0] = a_constant * x[0];
	}
};
Dynamics dyn;

int
main()
{
	real_t t;
	real_t x[x_dim];
	//* integration loop 
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, &t, x);

	return 0;
}