#include "rk4_solver.hpp"
#include <iostream>
#include <chrono>

constexpr uint_t t_dim = 1e9;
constexpr uint_t x_dim = 4;
constexpr real_t h = 1. / (t_dim - 1);

real_t t = 0;
real_t x[x_dim] = { 0 };
uint_t i;

//* dt__x = f(t, x) = t, x(0) = 0
//* x = a*t^2
void
ode_fun(const real_t t, const real_t[], const uint_t, real_t dt__x[])
{
	dt__x[0] = t;
	dt__x[1] = .5*t;
	dt__x[2] = 2*t;
	dt__x[3] = .25*t;
}

int
main()
{
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	rk4_solver::loop<ode_fun, t_dim, x_dim>(t, x, h, &t, x);

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << " time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
		  << std::endl;

	return 0;
}
