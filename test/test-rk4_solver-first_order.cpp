#include "matrix_rw.hpp"
#include "rk4_solver.hpp"
#include <cmath>

using Uint_T = rk4_solver::Uint_T;
using Real_T = rk4_solver::Real_T;

//* setup
const std::string dat_dir = "../dat";
const std::string test_name = "test-rk4_solver-first_order";
const std::string dat_prefix = dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string x_arr_fname = "x_arr.dat";

constexpr Uint_T sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 1.;
constexpr Uint_T t_dim = sample_freq*(t_final - t_init) + 1;
constexpr Uint_T x_dim = 1;
constexpr Real_T x_init[x_dim] = {1.};
constexpr Real_T a_const = 1.;

#ifdef __USE_SINGLE_PRECISION__
constexpr Real_T error_thres = 1e-5;
#else
constexpr Real_T error_thres = 1e-9;
#endif

struct Dynamics {
	//* dt__x = f(t, x) = a*x
	//* x = exp(a*t)
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const Uint_T, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = a_const * x[0];
	}
};
Dynamics dyn;

int
main()
{
	//* test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, &t, x);
	rk4_solver::cum_loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, t_arr, x_arr);

	//* write test data
	matrix_rw::write<t_dim, 1>(dat_prefix + t_arr_fname, t_arr);
	matrix_rw::write<t_dim, x_dim>(dat_prefix + x_arr_fname, x_arr);

	//* verify
	Real_T max_error = 0.;

	for (Uint_T i = 0; i < t_dim; ++i) {
		const Real_T(&x_)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr);
		const Real_T error = std::abs(x_[0] - std::exp(a_const * t_arr[i]));

		if (error > max_error) {
			max_error = error;
		}
	}

	//* loop vs cum_loop sanity check
	Real_T max_loop_error = 0.;
	const Real_T(&x_final)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(t_dim - 1, x_arr);

	for (Uint_T i = 0; i < x_dim; ++i) {
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
