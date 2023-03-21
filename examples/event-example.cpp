#include "rk4_solver/loop.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e5;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 2.;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 2;
constexpr Real_T x_init[x_dim] = {1., 0.};
constexpr Real_T e_restitution = .75;
constexpr Real_T gravity_const = 9.806;

struct Dynamics {
	/*
	 * Ball in vertical axis:
	 * dt_x =  [x2; -g]
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = x[1];
		dt_x[1] = -gravity_const;
	}
	/*
	 * Impact event:
	 * x+ = 0
	 * dt_x+ = -e * dt_x
	 */
	bool
	event_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&x_plus)[x_dim])
	{
		bool did_occur = false;

		if (x[0] <= 0) {
			x_plus[0] = 0;
			x_plus[1] = -e_restitution * x[1];
			did_occur = true;
		}
		return did_occur;
	}
};
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	rk4_solver::Event<x_dim, Dynamics> event(dynamics, &Dynamics::event_fun);

	Real_T t;
	Real_T x[x_dim];
	rk4_solver::loop<t_dim>(integrator, event, t_init, x_init, t, x, true);

	return 0;
}
