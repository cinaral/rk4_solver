#include "rk4_solver.hpp"
#include <chrono>
#include <iostream>

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e8;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 1.;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 4;

struct Dynamics {
	/*
	 * dt_x = f(t, x) = t, x(0) = 0
	 * x = v*t^2
	 */
	void
	ode_fun(const Real_T t, const Real_T (&)[x_dim], const size_t, Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = t;
		dt_x[1] = .5 * t;
		dt_x[2] = 2 * t;
		dt_x[3] = .25 * t;
	}
};
Dynamics dyn;

void print_elapsed_since(const std::chrono::time_point<std::chrono::high_resolution_clock> &start);

int
main()
{
	Real_T t;
	Real_T x[x_dim];
	const auto start = std::chrono::high_resolution_clock::now();
	rk4_solver::loop<t_dim>(dyn, &Dynamics::ode_fun, t_init, x, time_step, &t, x);
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