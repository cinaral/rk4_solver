#include "rk4_solver/loop.hpp"

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;


constexpr uint_t sample_freq = 1e5;
constexpr real_t time_step = 1. / sample_freq;
constexpr real_t t_init = 0.;
constexpr real_t t_final = 2.;
constexpr uint_t t_dim = sample_freq*(t_final - t_init) + 1;
constexpr uint_t x_dim = 2;
constexpr real_t x_init[x_dim] = {1., 0.};
constexpr real_t e_restitution = .75;
constexpr real_t gravity_const = 9.806;


struct Dynamics {
	//* Ball in vertical axis:
	//* dt__x =  [x2; -g]
	void
	ode_fun(const real_t, const real_t (&x)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
	{
		dt__x[0] = x[1];
		dt__x[1] = -gravity_const;
	}
	//* Impact event:
	//* x+ = 0
	//* dt__x+ = -e * dt__x
	bool
	event_fun(const real_t, const real_t (&x)[x_dim], const uint_t, real_t (&x_plus)[x_dim])
	{
		if (x[0] <= 0) {
			x_plus[0] = 0;
			x_plus[1] = -e_restitution * x[1];
		}
		//* don't stop the integration
		return false;
	}
};
Dynamics dyn;

int
main()
{
	real_t t;
	real_t x[x_dim];
	//* integration loop with events
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t_init, x_init, time_step, &t, x);

	return 0;
}
