#include "test_config.hpp"

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

//* setup
const std::string test_name = "sine-test";
const std::string dat_prefix = test_config::dat_dir + "/" + test_name + "-";

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0;
constexpr Real_T t_final = 1;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x_init[x_dim] = {0};
constexpr Real_T sine_freq = 5.;

#ifdef __USE_SINGLE_PRECISION__
constexpr Real_T error_thres = 1e-5;
#else
constexpr Real_T error_thres = 1e-9;
#endif

struct Dynamics {
	/*
	 * dt__x = f(t, x) = 2*pi*f*cos(t*2*pi*f)
	 * x = sin(t*2*pi*f)
	 */
	void
	ode_fun(const Real_T t, const Real_T (&)[x_dim], const size_t, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = 2 * M_PI * sine_freq * cos(t * 2 * M_PI * sine_freq);
	}
};
Dynamics dyn;

int
main()
{
	//* 1. read the reference data
	//* no reference data

	//* 2. test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	rk4_solver::loop<t_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, &t, x);
	rk4_solver::cum_loop<t_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, t_arr,
	                            x_arr);

	//* 3. write the test data
	matrix_rw::write<t_dim, 1>(dat_prefix + test_config::t_arr_fname, t_arr);
	matrix_rw::write<t_dim, x_dim>(dat_prefix + test_config::x_arr_fname, x_arr);

	//* 4. verify the results
	Real_T x_arr_chk[t_dim * x_dim];

	for (size_t i = 0; i < t_dim; ++i) {
		Real_T x_chk[x_dim] = {std::sin(t_arr[i] * 2 * M_PI * sine_freq)};
		matrix_op::replace_row<t_dim, x_dim>(i, x_chk, x_arr_chk);
	}
	Real_T max_error = test_config::compute_max_error<t_dim, x_dim>(x_arr, x_arr_chk);

	//* loop vs cum_loop sanity check
	Real_T max_loop_error = 0.;
	const Real_T(&x_final)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(t_dim - 1, x_arr);

	for (size_t i = 0; i < x_dim; ++i) {
		const Real_T error = std::abs(x_final[i] - x[i]);

		if (error > max_loop_error) {
			max_loop_error = error;
		}
	}

	if (max_error < error_thres && max_loop_error < std::numeric_limits<Real_T>::epsilon()) {
		return 0;
	} else {
		return 1;
	}
}
