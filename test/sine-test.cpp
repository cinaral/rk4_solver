#include "test_config.hpp"

//* setup
const std::string test_name = "sine-test";
const std::string dat_prefix = test_config::dat_dir + "/" + test_name + "-";

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 1.;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 1;
constexpr Real_T x_init[x_dim] = {0.};
constexpr Real_T sine_freq = 5.;

#ifdef USE_SINGLE_PRECISION
constexpr Real_T error_thres = 1e-5;
#else
constexpr Real_T error_thres = 1e-9;
#endif

struct Dynamics {
	/*
	 * dt_x = f(t, x) = 2*pi*f*cos(t*2*pi*f)
	 * x = sin(t*2*pi*f)
	 */
	void
	ode_fun(const Real_T t, const Real_T (&)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = 2 * M_PI * sine_freq * cos(t * 2 * M_PI * sine_freq);
	}
};
Dynamics dynamics;

int
main()
{
	//* 1. read the reference data
	//* no reference data

	//* 2. test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[1][t_dim];
	Real_T x_arr[t_dim][x_dim];

	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	rk4_solver::loop<t_dim>(integrator, t_init, x_init, t, x);

	integrator.reset();
	rk4_solver::loop(integrator, t_init, x_init, t_arr[0], x_arr);

	//* 3. write the test data
	matrix_rw::write(dat_prefix + test_config::t_arr_fname, t_arr);
	matrix_rw::write(dat_prefix + test_config::x_arr_fname, x_arr);

	//* 4. verify the results
	Real_T x_arr_ref[t_dim][x_dim];

	for (size_t i = 0; i < t_dim; ++i) {
		Real_T x_ref[x_dim] = {
		    static_cast<Real_T>(std::sin(t_arr[0][i] * 2 * M_PI * sine_freq))};
		matrix_op::replace_row(i, x_ref, x_arr_ref);
	}
	Real_T max_error = test_config::compute_max_error(x_arr, x_arr_ref);

	//* cumulative loop sanity check
	Real_T max_loop_error = 0.;
	const Real_T(&x_final)[x_dim] = x_arr[t_dim - 1];

	for (size_t i = 0; i < x_dim; ++i) {
		const Real_T error = std::abs(x_final[i] - x[i]);

		if (error > max_loop_error) {
			max_loop_error = error;
		}
	}

	if (max_error < error_thres && max_loop_error <= std::numeric_limits<Real_T>::epsilon()) {
		return 0;
	} else {
		printf("max_error = %.3g\n", max_error);
		printf("max_loop_error = %.3g\n", max_loop_error);
		return 1;
	}
}
