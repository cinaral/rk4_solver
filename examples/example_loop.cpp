#include "rk4_solver.hpp"
#include "types.hpp"

constexpr uint_t t_dim = 1e3;
constexpr uint_t x_dim = 1;

//* dt__x = f(t, x) = t, x(0) = 0
//* x = 1/2*t^2
void
ode_fun(const real_t t, const real_t[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = t;
}

real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];

int
main()
{
	//* initialize t_arr
	for (uint_t i = 0; i < t_dim; i++) {
		t_arr[i] =  1. * i / (t_dim - 1);
	}

	//* integration loop
	rk4_solver::loop<ode_fun, t_dim, x_dim>(t_arr, x_arr);

	return 0;
}