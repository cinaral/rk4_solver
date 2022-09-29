#include "rk4_solver/step.hpp"

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;

constexpr uint_t x_dim = 2;
constexpr real_t u = 1e3;
constexpr real_t h = 1e-3;
constexpr real_t t = 0;
constexpr real_t x[x_dim] = {0, 0};
constexpr uint_t i = 0;

struct Dynamics {
	//* dt__x = f(t, x) = [x[1]; u], x(0) = [0; 0];
	void
	ode_fun(const real_t, const real_t (&x)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
	{
		dt__x[0] = x[1];
		dt__x[1] = u;
	}
};
Dynamics dyn;

int
main()
{
	real_t x_next[x_dim];
	//* integration step
	rk4_solver::step<Dynamics, x_dim>(dyn, &Dynamics::ode_fun, t, x, h, i, x_next);

	return 0;
}