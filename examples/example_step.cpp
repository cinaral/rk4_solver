#include "rk4_solver.hpp"
#include "types.hpp"

constexpr uint_t x_dim = 2;
constexpr real_t u = 1e3;

//* dt__x = f(t, x) = [x[1]; u], x(0) = [0; 0];
void
ode_fun(const real_t, const real_t x[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = x[1];
	OUT_dt__x[1] = u;
}

real_t x[x_dim] = {0, 0};
real_t h = 1e-3;
real_t OUT_x_next[x_dim];

int
main()
{
	//* integration step
	rk4_solver::step<ode_fun, x_dim>(0, x, 0, h, OUT_x_next);

	return 0;
}