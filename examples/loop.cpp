#include "rk4_solver/cum_loop.hpp"

using size_t = rk4_solver::size_t;
using Real_T = rk4_solver::Real_T;

constexpr size_t sample_freq = 1e3;
constexpr Real_T h = 1. / sample_freq;
constexpr Real_T t0 = 0.;
constexpr Real_T tf = 1.;
constexpr size_t t_dim = sample_freq*(tf - t0) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x0[x_dim] = {1.};

Real_T t_arr[t_dim];
Real_T x_arr[t_dim * x_dim];

struct Dynamics {
	//* dt__x = f(t, x) = t, x(0) = 0
	//* x = 1/2*t^2
	void
	ode_fun(const Real_T t, const Real_T (&)[x_dim], const size_t, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = t;
	}
};
Dynamics dyn;

int
main()
{
	//* integration loop with cumulatively saved data arrays
	rk4_solver::cum_loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t0, x0, h, t_arr, x_arr);

	return 0;
}
