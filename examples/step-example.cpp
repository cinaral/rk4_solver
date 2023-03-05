#include "rk4_solver/step.hpp"

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t x_dim = 2;
constexpr Real_T u = 1e3;
constexpr Real_T h = 1e-3;
constexpr Real_T t = 0;
constexpr Real_T x[x_dim] = {0, 0};
constexpr size_t i = 0;

struct Dynamics {
	/*
	 * dt_x = f(t, x) = [x[1]; u], x(0) = [0; 0];
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = x[1];
		dt_x[1] = u;
	}
};
Dynamics dyn;

int
main()
{
	Real_T x_next[x_dim];
	//* integration step
	rk4_solver::step(dyn, &Dynamics::ode_fun, t, x, h, i, x_next);
	
	return 0;
}
