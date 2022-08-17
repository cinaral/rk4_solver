#include "rk4_solver.hpp"
#include "types.hpp"

constexpr uint_t t_dim = 1e3;
constexpr uint_t x_dim = 1;
constexpr real_t t0 = 0.;
constexpr real_t tf = 1.;
constexpr real_t x0[x_dim] = {1.};
constexpr real_t h = tf / (t_dim - 1);

real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];

//* dt__x = f(t, x) = t, x(0) = 0
//* x = 1/2*t^2
void
ode_fun(const real_t t, const real_t[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = t;
}

int
main()
{
	//* integration loop
	rk4_solver::cum_loop<ode_fun, t_dim, x_dim>(t0, x0, h, t_arr, x_arr);

	return 0;
}