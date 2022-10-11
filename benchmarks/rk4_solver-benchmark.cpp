#include "rk4_solver.hpp"
#include <chrono>
#include <iostream>

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;

constexpr uint_t sample_freq = 1e9;
constexpr real_t time_step = 1. / sample_freq;
constexpr real_t t_init = 0.;
constexpr real_t t_final = 1.;
constexpr uint_t t_dim = sample_freq*(t_final - t_init) + 1;
constexpr uint_t x_dim = 4;

struct Dynamics {
	//* dt__x = f(t, x) = t, x(0) = 0
	//* x = v*t^2
	void
	ode_fun(const real_t t, const real_t (&)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
	{
		dt__x[0] = t;
		dt__x[1] = .5 * t;
		dt__x[2] = 2 * t;
		dt__x[3] = .25 * t;
	}
};
Dynamics dyn;

void
print_elapsed_since(const std::chrono::time_point<std::chrono::high_resolution_clock> &start);

int
main()
{
	real_t t;
	real_t x[x_dim];
	const auto start = std::chrono::high_resolution_clock::now();
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x, time_step, &t, x);
	print_elapsed_since(start);

	return 0;
}

void
print_elapsed_since(const std::chrono::time_point<std::chrono::high_resolution_clock> &start)
{
	const auto end = std::chrono::high_resolution_clock::now();
	const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "elapsed (ms): " << diff.count() << std::endl;
}