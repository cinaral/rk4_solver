#include "rk4_solver/integrator.hpp"
#include <chrono>
#include <cstdio>

using rk4_solver::Real_T;
using rk4_solver::size_t;

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr size_t x_dim = 3;
constexpr Real_T t_init = 0.;
constexpr Real_T x_init[x_dim] = {1e1, 1e0, 0};
constexpr size_t benchmark_duration_ms = 5;

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
	Real_T t = t_init;
	Real_T x[x_dim];
	matrix_op::replace_row<1>(0, x_init, x); //* initialize x

	printf("Integrating 3rd order linear ODE for %zu ms... ", benchmark_duration_ms);
	auto sample_tp = std::chrono::high_resolution_clock::now();
	auto now_tp = std::chrono::high_resolution_clock::now();
	auto since_sample = sample_tp - now_tp;

	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);

	while (true) {
		integrator.step(t, x, t, x);

		now_tp = std::chrono::high_resolution_clock::now();
		since_sample = now_tp - sample_tp;

		if (since_sample >= std::chrono::milliseconds(benchmark_duration_ms)) {
			const size_t step_count = integrator.get_step_count();
			printf("Done.\nx at t = %.3g s: [%.3g; %.3g; %.3g]\n", t, x[0], x[1], x[2]);
			printf("Score: %.3g steps per second (%.3g steps)\n",
			       static_cast<Real_T>(step_count) / benchmark_duration_ms * 1e3,
			       static_cast<Real_T>(step_count));
			sample_tp = std::chrono::high_resolution_clock::now();
			break;
		}
	}
	return 0;
}
