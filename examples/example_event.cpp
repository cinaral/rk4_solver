#include "rk4_solver.hpp"

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;

constexpr uint_t t_dim = 1e5;
constexpr uint_t x_dim = 2;
constexpr real_t t0 = 0.;
constexpr real_t tf = 2.;
constexpr real_t x0[x_dim] = {1., 0.};
constexpr real_t h = tf / (t_dim - 1);

//* Ball equations:
//* dt__x =  [x2; -g]
constexpr real_t e = .75;
constexpr real_t g = 9.806;

real_t t = 0;
real_t x[x_dim];

void
ode_fun(const real_t, const real_t (&x)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
{
	dt__x[0] = x[1];
	dt__x[1] = -g;
}

bool
event_fun(const real_t, const real_t (&x)[x_dim], const uint_t, real_t (&x_plus)[x_dim])
{
	if (x[0] <= 0) {
		x_plus[0] = 0;
		x_plus[1] = -e * x[1];
	}

	//* don't stop the integration
	return false;
}

int
main()
{
	//* integration loop with events
	rk4_solver::loop<t_dim, x_dim, ode_fun, event_fun>(t0, x0, h, &t, x);

	return 0;
}
