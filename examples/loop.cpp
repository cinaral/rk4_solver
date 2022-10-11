#include "rk4_solver/cum_loop.hpp"

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;

constexpr uint_t sample_freq = 1e3;
constexpr real_t time_step = 1. / sample_freq;
constexpr real_t t_init = 0.;
constexpr real_t t_final = 1.;
constexpr uint_t t_dim = sample_freq*(t_final - t_init) + 1;
constexpr uint_t x_dim = 1;
constexpr real_t x_init[x_dim] = {1.};

real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];

struct Dynamics {
	//* dt__x = f(t, x) = t, x(0) = 0
	//* x = 1/2*t^2
	void
	ode_fun(const real_t t, const real_t (&)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
	{
		dt__x[0] = t;
	}
};
Dynamics dyn;

int
main()
{
	//* integration loop with cumulatively saved data arrays
	rk4_solver::cum_loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, t_arr, x_arr);

	return 0;
}