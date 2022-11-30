#include "rk4_solver/loop.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e5;
constexpr Real_T h = 1. / sample_freq;
constexpr Real_T t0 = 0.;
constexpr Real_T tf = 2.;
constexpr size_t t_dim = sample_freq * (tf - t0) + 1;
constexpr size_t x_dim = 2;
constexpr Real_T x0[x_dim] = {1., 0.};
constexpr Real_T e_restitution = .75;
constexpr Real_T gravity_const = 9.806;

struct Dynamics {
	//* Ball in vertical axis:
	//* dt__x =  [x2; -g]
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = x[1];
		dt__x[1] = -gravity_const;
	}
	//* Impact event:
	//* x+ = 0
	//* dt__x+ = -e * dt__x
	bool
	event_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&x_plus)[x_dim])
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
	Real_T t;
	Real_T x[x_dim];
	//* integration loop with events
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t0,
	                                         x0, h, &t, x);

	return 0;
}
