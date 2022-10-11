#include "matrix_rw.hpp"
#include "rk4_solver.hpp"
#include <cmath>

using uint_t = rk4_solver::uint_t;
using real_t = rk4_solver::real_t;

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

//* setup
const std::string dat_dir = "../dat";
const std::string test_name = "test-rk4_solver-sine";
const std::string dat_prefix = dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string x_arr_fname = "x_arr.dat";

constexpr uint_t sample_freq = 1e3;
constexpr real_t time_step = 1. / sample_freq;
constexpr real_t t_init = 0;
constexpr real_t t_final = 1;
constexpr uint_t t_dim = sample_freq*(t_final - t_init) + 1;
constexpr uint_t x_dim = 1;
constexpr real_t x_init[x_dim] = {0};
constexpr real_t sine_freq = 5.;

#ifdef __USE_SINGLE_PRECISION__
constexpr real_t error_thres = 1e-5;
#else
constexpr real_t error_thres = 1e-9;
#endif

struct Dynamics {
	//* dt__x = f(t, x) = 2*pi*f*cos(t*2*pi*f)
	//* x = sin(t*2*pi*f)
	void
	ode_fun(const real_t t, const real_t (&)[x_dim], const uint_t, real_t (&dt__x)[x_dim])
	{
		dt__x[0] = 2 * M_PI * sine_freq * cos(t * 2 * M_PI * sine_freq);
	}
};
Dynamics dyn;

int
main()
{
	//* test
	real_t t = 0;
	real_t x[x_dim];
	real_t t_arr[t_dim];
	real_t x_arr[t_dim * x_dim];
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, &t, x);
	rk4_solver::cum_loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, t_arr, x_arr);

	//* write test data
	matrix_rw::write<t_dim, 1>(dat_prefix + t_arr_fname, t_arr);
	matrix_rw::write<t_dim, x_dim>(dat_prefix + x_arr_fname, x_arr);

	//* verify
	real_t max_error = 0.;

	for (uint_t i = 0; i < t_dim; ++i) {
		const real_t(&x_)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr);
		const real_t error = std::abs(x_[0] - std::sin(t_arr[i] * 2 * M_PI * sine_freq));

		if (error > max_error) {
			max_error = error;
		}
	}

	//* loop vs cum_loop sanity check
	real_t max_loop_error = 0.;
	const real_t(&x_final)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(t_dim - 1, x_arr);

	for (uint_t i = 0; i < x_dim; ++i) {
		const real_t error = std::abs(x_final[i] - x[i]);

		if (error > max_loop_error) {
			max_loop_error = error;
		}
	}

	if (max_error < error_thres && max_loop_error < std::numeric_limits<real_t>::epsilon()) {
		return 0;
	} else {
		return 1;
	}
}
