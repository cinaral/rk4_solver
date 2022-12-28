#include "test_config.hpp"

/*
 * See https://www.mathworks.com/help/matlab/math/ode-event-location.html
 * for the official MATLAB example counterpart of what we are trying to achieve in this test.
 */

//* setup
const std::string test_name = "ball-test";
const std::string dat_prefix = test_config::dat_dir + "/" + test_name + "-";
const std::string ref_dat_prefix = test_config::ref_dat_dir + "/" + test_name + "-";

constexpr size_t sample_freq = 1e4;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 2.;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 2;
constexpr Real_T x_init[x_dim] = {1., 0.};
constexpr Real_T e_restitution = .75;
constexpr Real_T gravity_const = 9.806;

#ifdef __USE_SINGLE_PRECISION__
constexpr Real_T error_thres = 2e-3;
#else
constexpr Real_T error_thres = 1e-3;
#endif
constexpr size_t verify_idx = 0;

struct Dynamics {
	/*
	 * Ball equations:
	 * dt__x =  [x2; -g]
	 */
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = x[1];
		dt__x[1] = -gravity_const;
	}

	bool
	event_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&x_plus)[x_dim])
	{
		if (x[0] <= 0) {
			x_plus[0] = 0;
			x_plus[1] = -e_restitution * x[1];
		}
		return false;
	}
};
Dynamics dyn;

int
main()
{
	//* 1. read the reference data
	Real_T x_arr_chk[t_dim * x_dim];
	matrix_rw::read<t_dim, x_dim>(ref_dat_prefix + test_config::x_arr_chk_fname, x_arr_chk);

	//* 2. test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	rk4_solver::loop<t_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t_init, x_init,
	                        time_step, &t, x);
	rk4_solver::cum_loop<t_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t_init, x_init,
	                            time_step, t_arr, x_arr);

	//* 3. write the test data
	matrix_rw::write<t_dim, 1>(dat_prefix + test_config::t_arr_fname, t_arr);
	matrix_rw::write<t_dim, x_dim>(dat_prefix + test_config::x_arr_fname, x_arr);

	//* 4. verify the results
	Real_T x_0_arr[t_dim * 1];
	Real_T x_0_arr_chk[t_dim * 1];

	for (size_t i = 0; i < t_dim; ++i) {
		const Real_T(&x)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr);
		const Real_T(&x_chk)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr_chk);
		x_0_arr[i] = x[verify_idx];
		x_0_arr_chk[i] = x_chk[verify_idx];
	}
	Real_T max_error = test_config::compute_max_error<t_dim, 1>(x_0_arr, x_0_arr_chk);

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
