#include "rk4_solver/loop.hpp"
#include <chrono>
#include <iostream>

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr size_t x_dim = 3;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 2e1;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr Real_T x_init[x_dim] = {1e1, 1e0, 0};

struct Dynamics {
	/*
	 * dt_x = A * x
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = x[1];
		dt_x[1] = x[2];
		dt_x[2] = -a0 * dt_x[0] - a1 * x[1] - a2 * x[2];
	}
	const Real_T a0 = 1e-1;
	const Real_T a1 = 1e-2;
	const Real_T a2 = 1e-3;
};
Dynamics dynamics;

int
main()
{
#ifdef DO_NOT_USE_HEAP
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
#else
	Real_T(&t_arr)[t_dim] = *(Real_T(*)[t_dim]) new Real_T[t_dim];
	Real_T(&x_arr)[t_dim][x_dim] = *(Real_T(*)[t_dim][x_dim]) new Real_T[t_dim][x_dim];
#endif

	printf("Cumulatively integrating 3rd order linear ODE for %.3g steps... ",
	       static_cast<Real_T>(t_dim));

	const auto start_tp = std::chrono::high_resolution_clock::now();
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);

	rk4_solver::loop<t_dim>(integrator, t_init, x_init, t_arr, x_arr);

	const auto now_tp = std::chrono::high_resolution_clock::now();
	const auto since_sample_ns =
	    std::chrono::duration_cast<std::chrono::nanoseconds>(now_tp - start_tp);

	const Real_T(&x_final)[x_dim] = x_arr[t_dim - 1];

	printf("Done.\nx at t = %.3g s: [%.3g; %.3g; %.3g]\n", t_arr[t_dim - 1], x_final[0],
	       x_final[1], x_final[2]);
	printf("Score: %.3g steps per second (%g ms)\n",
	       static_cast<Real_T>(t_dim) / since_sample_ns.count() * 1e9,
	       static_cast<Real_T>(since_sample_ns.count()) / 1e9);

	return 0;
}
